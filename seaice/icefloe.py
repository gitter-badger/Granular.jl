#!/usr/bin/env python
import numpy


class IceFloe:
    def init(self,
             lin_pos=[0., 0.],
             lin_vel=[0., 0.],
             lin_acc=[0., 0.],
             force=[0., 0.],
             ang_pos=0.,
             ang_vel=0.,
             ang_acc=0.,
             torque=0.,
             radius=1.,
             thickness=1.,
             density=934.):
        '''
        Initializing function with modifiable default values
        :param lin_pos: Floe linear position [m]
        :type lin_pos: list or numpy.array
        :param lin_vel: Floe linear velocity [m/s]
        :type lin_vel: list or numpy.array
        :param lin_acc: Floe linear acceleration [m/s^2]
        :type lin_acc: list or numpy.array
        :param force: Sum of forces [N]
        :type force: list or numpy.array
        :param ang_pos: Floe angular position [rad]
        :type ang_pos: float
        :param ang_vel: Floe angular velocity [rad/s]
        :type ang_vel: float
        :param ang_acc: Floe angular acceleration [rad/s^2]
        :type ang_acc: float
        :param torque: Sum of forces [N]
        :type torque: float
        :param radius: Floe radius [m]
        :type radius: float
        :param thickness: Floe thickness [m]
        :type thickness: float
        :param density: Floe density [kg/m^3]
        :type density: float
        '''

        self.lin_pos = numpy.array(lin_pos)
        self.lin_vel = numpy.array(lin_vel)
        self.lin_acc = numpy.array(lin_acc)
        self.force = numpy.array(force)
        self.ang_pos = numpy.array(ang_pos)
        self.ang_vel = numpy.array(ang_vel)
        self.ang_acc = numpy.array(ang_acc)
        self.force = numpy.array(force)
        self.radius = radius
        self.thickness = thickness
        self.density = density

    def update_position(dt):
        pass
