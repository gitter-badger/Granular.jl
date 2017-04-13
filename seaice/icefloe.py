#!/usr/bin/env python
import numpy


class IceFloeCylindrical:
    '''
    Cylindrical ice floe object.

    :param lin_pos: Floe linear position [m]
    :type lin_pos: list or numpy.array
    :param thickness: Floe thickness [m]
    :type thickness: float
    :param contact_radius: Floe radius during interactions [m]
    :type contact_radius: float
    :param areal_radius: Floe areal radius on the sea surface [m].  If not
    set, this parameter will equal the `contact_radius` value.
    :type areal_radius: float
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
    :param density: Floe density [kg/m^3]
    :type density: float
    :param rotating: The floe is free to rotate
    :type rotating: bool
    :param fixed: The floe is free to move linearly and/or rotationally
    :type fixed: bool
    '''
    def __init__(self,
             lin_pos,
             thickness,
             contact_radius,
             areal_radius=None,
             lin_vel=[0., 0.],
             lin_acc=[0., 0.],
             force=[0., 0.],
             ang_pos=0.,
             ang_vel=0.,
             ang_acc=0.,
             torque=0.,
             density=934.,
             rotating=True,
             fixed=False):

        self.lin_pos = numpy.array(lin_pos)
        self.thickness = thickness
        self.contact_radius = contact_radius
        if areal_radius is None:
            self.areal_radius = contact_radius
        self.lin_vel = numpy.array(lin_vel)
        self.lin_acc = numpy.array(lin_acc)
        self.force = numpy.array(force)
        self.ang_pos = numpy.array(ang_pos)
        self.ang_vel = numpy.array(ang_vel)
        self.ang_acc = numpy.array(ang_acc)
        self.force = numpy.array(force)
        self.density = density
        self.rotating = rotating
        self.fixed = fixed

    def surface_area(self):
        '''
        Determine the current floe surface area.

        :returns: The floe mass based on its areal radius.
        :return type: float
        '''
        return numpy.pi*self.areal_radius**2.

    def volume(self):
        '''
        Determine the current floe volume.

        :returns: The floe mass based on its areal radius, thickness.
        :return type: float
        '''
        return self.thickness*self.surface_area()

    def mass(self):
        '''
        Determine the current floe mass.

        :returns: The floe mass based on its areal radius, thickness, and
        density.
        :return type: float
        '''
        return self.density*self.volume()

    def moment_of_inertia_vertical(self):
        '''
        Determines the rotational moment of inertia for rotation around the
        floe center with a vertical rotation axis.

        :returns: The vertical rotational moment of inertia
        :return type: float
        '''
        return 0.5*self.mass()*self.areal_radius**2.

    def find_accelerations(self):
        '''
        Determine current linear and angular accelerations based on the current
        sum of forces.
        '''
        if self.fixed is False:
            self.lin_acc = self.force/self.mass()
            if self.rotating:
                self.ang_acc = self.torque/self.moment_of_inertia_vertical()

    def update_position(self, dt, method='TY3'):
        '''
        Update the kinematics through explicit temporal integration using a
        third-order Taylor expansion.

        :param dt: The time step length
        :type dt: float
        :param method: The integration method to choose ('TY2' or 'TY3'
            (default))
        :type str: str
        '''
        if method == 'TY2':
            self.update_position_TY2(dt)

        elif method == 'TY3':
            self.update_position_TY3(dt)

        else:
            raise Exception('Error: Intergration method not understood.')

    def update_position_TY2(self, dt):
        '''
        Update the kinematics through explicit temporal integration of Newton's
        second law using a second-order Taylor expansion.

        :param dt: The time step length
        :type dt: float
        '''

        self.find_accelerations()

        self.lin_pos += self.lin_vel*dt + 0.5*self.lin_acc*dt**2.
        self.ang_pos += self.ang_vel*dt + 0.5*self.ang_acc*dt**2.

        self.lin_vel += self.lin_acc*dt
        self.ang_vel += self.ang_acc*dt

    def update_position_TY3(self, dt):
        '''
        Update the kinematics through explicit temporal integration of Newton's
        second law using a third-order Taylor expansion.

        :param dt: The time step length
        :type dt: float
        '''
        lin_acc0 = self.lin_acc
        ang_acc0 = self.ang_acc
        self.find_accelerations()

        d_lin_acc_dt = (self.lin_acc - lin_acc0)/dt
        d_ang_acc_dt = (self.ang_acc - ang_acc0)/dt

        self.lin_pos += self.lin_vel*dt + 0.5*self.lin_acc*dt**2. \
            + 1./6.*d_lin_acc_dt*dt**3.
        self.ang_pos += self.ang_vel*dt + 0.5*self.ang_acc*dt**2. \
            + 1./6.*d_ang_acc_dt*dt**3.

        self.lin_vel += self.lin_acc*dt + 0.5*d_lin_acc_dt*dt**2.
        self.ang_vel += self.ang_acc*dt + 0.5*d_ang_acc_dt*dt**2.
