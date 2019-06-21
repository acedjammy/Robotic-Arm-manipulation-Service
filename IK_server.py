#!/usr/bin/env python

# Copyright (C) 2017 Udacity Inc.
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
# All Rights Reserved.
# Author: Yang Ding

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from sympy import *

#Create the subfunction to generate homogeneous transformation based on DH parameters
def TF_Mat(alpha, a, d, q):
    TF = Matrix([[            cos(q),           -sin(q),           0,             a],
                 [ sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
                 [ sin(q)*sin(alpha), cos(q)*sin(alpha),  cos(alpha),  cos(alpha)*d],
                 [                 0,                 0,           0,             1]])
    return TF

# IK handler service
def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')                                 # joint angles theta
        d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')                                 # link offsets
        a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')                                 # link lengths
        alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7') # joint twist angles

        dh = {alpha0:      0, a0:      0, d1:  0.75, q1:        q1,
              alpha1: -pi/2., a1:   0.35, d2:     0, q2: -pi/2.+q2,
              alpha2:      0, a2:   1.25, d3:     0, q3:        q3,
              alpha3: -pi/2., a3: -0.054, d4:   1.5, q4:        q4,
              alpha4:  pi/2., a4:      0, d5:     0, q5:        q5,
              alpha5: -pi/2., a5:      0, d6:     0, q6:        q6,
              alpha6:      0, a6:      0, d7: 0.303, q7:         0}

        # Substitute DH_Table
        T01 = TF_Mat(alpha0, a0, d1, q1).subs(dh)
        T12 = TF_Mat(alpha1, a1, d2, q2).subs(dh)
        T23 = TF_Mat(alpha2, a2, d3, q3).subs(dh)
        T34 = TF_Mat(alpha3, a3, d4, q4).subs(dh)
        T45 = TF_Mat(alpha4, a4, d5, q5).subs(dh)
        T56 = TF_Mat(alpha5, a5, d6, q6).subs(dh)
        T67 = TF_Mat(alpha6, a6, d7, q7).subs(dh)

        # Transform from base link to link 3 and end effector
        T03 = T01 * T12 * T23
        T07 = T01 * T12 * T23 * T34 * T45 * T56 * T67 

        
        # Find EE rotation matrix RPY (Roll, Pitch, Yaw)
        r,p,y = symbols('r p y')

        # Roll
        Rx = Matrix([[       1,       0,       0],
                     [       0,  cos(r), -sin(r)],
                     [       0,  sin(r),  cos(r)]])
        # Pitch
        Ry = Matrix([[  cos(p),       0,  sin(p)],
                     [       0,       1,       0],
                     [ -sin(p),       0,  cos(p)]])
        # Yaw
        Rz = Matrix([[  cos(y), -sin(y),       0],
                     [  sin(y),  cos(y),       0],
                     [       0,       0,       1]])

        #Rotation matrix of the Euler Angles represented in symbols
        REE = Rz * Ry * Rx
        
        # Compensate for rotation discrepancy between DH parameters and Gazebo
        # Correction for the rotation
        Rcorr = Matrix([[0, 0, 1],[0, -1, 0],[1, 0, 0]])
        REE = REE * Rcorr     
        
        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

            # Extract end-effector position and orientation from request
            # px,py,pz = end-effector position
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z
            
            # store EE position in a matrix
            EE = Matrix([[px],
                        [py],
                        [pz]])
            
            # roll, pitch, yaw from quaternion
            (roll,pitch,yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x,
                 req.poses[x].orientation.y,
                 req.poses[x].orientation.z,
                 req.poses[x].orientation.w])

            #Substitute value
            REE = REE.subs({'r': roll, 'p': pitch, 'y': yaw})

            # Calculate Wrist Center location with respect to base link
            WC = EE - (0.303) * REE[:,2]

            # Calculate joint angles using Geometric IK method
            # Calculate theta1
            theta1 = atan2(WC[1],WC[0])
   
            #Calculate the sides and angles of the triangle
	    a = 1.25;
	    b = sqrt(1.5*1.5 + 0.054*0.054);
	    c = sqrt(pow((sqrt(WC[0]*WC[0] + WC[1]*WC[1]) - 0.35), 2) + pow((WC[2] - 0.75), 2))
	    A = acos((a*a - b*b -c*c)/(-2*b*c));
	    B = acos((b*b - a*a -c*c)/(-2*a*c));
	    C = acos((c*c - a*a -b*b)/(-2*a*b));

            # Calculate theta2 and theta3
	    theta2 = pi/2 - B - atan2(WC[2]-0.75, sqrt(WC[0]*WC[0]+WC[1]*WC[1])-0.35)
	    theta3 = pi/2 - C - atan2(0.054,1.5);

            #Calculate the roation matrix between frame 3 and end effector
            R03 = T03[0:3,0:3]
            R03 = R03.evalf(subs={q1: theta1, q2: theta2, q3:theta3})
            R36 = R03.transpose() * REE

            # Calculate theta4, theta5 and theta6
            theta5 = atan2(sqrt(R36[0,2]*R36[0,2] + R36[2,2]*R36[2,2]),R36[1,2])           
            theta4 = atan2(R36[2,2], -R36[0,2])
            theta6 = atan2(-R36[1,1],R36[1,0])

            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
            joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
            joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
