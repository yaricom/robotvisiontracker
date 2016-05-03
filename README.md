## Overview
Firm A is building a next generation robotics platform that will change the game in field service operations including asset 
inspection and repair. Firm A has defined a host of high value use cases and applications across industry that will support 
field engineers and other industrial workers be more productive and, more importantly, perform their jobs safely.

For one example high value use case, the company would like for a robot to detect and track a freight railcar brake release 
handle, the object of interest (OOI), so that the robot can grasp the handle.

The task is to develop an algorithm that can detect and track the OOI in video frames. The OOI is typically made of 0.5 inch 
round steel rod, bent to form a handle.

It was done as part of crowdsourcing contests on TopCoder: [Contest: Robot Vision Tracker](https://community.topcoder.com/longcontest/stats/?module=ViewOverview&rd=16405)
and [Contest: Robot Vision Tracker Extended ](https://community.topcoder.com/longcontest/stats/?module=ViewOverview&rd=16430)

### Special conditions ###
Additional requirement was to create crafted algorithm which may be compilled and runned in the AWS runner with specified 
limitation:
-	Time limit is 60 minutes per test case for training and 3 minutes for testing and the memory limit is 4096MB.
-	There is no explicit code size limit. The implicit source code size limit is around 1 MB (it is not advisable to submit codes of size close to that or larger).
-	The compilation time limit is 60 seconds. You can find information about compilers that we use, compilation options and processing server specifications here.

The main algorithm runner implemented in: [RobotVisionTracker.h](https://github.com/yaricom/robotvisiontracker/blob/master/RobotVisionTracker/RobotVisionTracker/RobotVisionTracker.h)
and [RobotVisionTrackerX.h](https://github.com/yaricom/robotvisiontracker/blob/master/RobotVisionTracker/RobotVisionTracker/RobotVisionTrackerX.h)
