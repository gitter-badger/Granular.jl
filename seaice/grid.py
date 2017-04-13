#!/usr/bin/env python
import numpy
import seaice


class SquareGrid:
    '''
    A two-dimensional, orthogonal and Cartesian grid, with options to add field
    values at the center, edges, or corners.

    :param nx: The number of grid cells along x.
    :type nx: int
    :param ny: The number of grid cells along y.
    :type ny: int
    '''
    def __init__(self, nx, ny, dx=None, dy=None, Lx=None, Ly=None):
        self.nx = int(nx)
        self.ny = int(ny)
        self.floes = []

    def setSize(self, origo=[0., 0.], dx=None, dy=None, Lx=None, Ly=None):
        '''
        Used to determine the spatial dimensions of the grid.  The user may
        provide cell widths (`dx` and/or `dy`) or grid lengths (`Lx` and/or
        `Ly`).  It is implied that the grid width equals the cell width
        multiplied by the number of cells along each dimension.

        :param origo: Shift the grid the value of this 2d vector.
        :type origo: numpy.array
        :param dx: Cell width along `x`.
        :type dx: float
        :param dy: Cell width along `y`.
        :type dy: float
        :param Lx: Grid width along `x`.
        :type dx: float
        :param Ly: Grid width along `y`.
        :type dy: float
        '''

        if dx is None and Lx is None:
            raise Exception('''Error: A SquareGrid must be initialized with cell
                            spacing or grid length along x.''')
        if dy is None and Ly is None:
            raise Exception('''Error: A SquareGrid must be initialized with cell
                            spacing or grid length along y.''')

        if dx:
            self.dx = float(dx)
            self.Lx = float(dx) * self.nx
        else:
            self.Lx = float(Lx)
            self.dx = float(Lx)/float(self.nx)
        if dy:
            self.dy = float(dy)
            self.Ly = float(dy) * self.ny
        else:
            self.Ly = float(Ly)
            self.dy = float(Ly)/float(self.ny)

        self.origo = numpy.array(origo)

    def getCenterCoordinate(self, i, j):
        '''
        Returns the center coordinate for the center of the cell with index `i`
        and `j`, sometimes referred to as the h-point.

        :param i: Cell index along `x`.
        :type i: int
        :param j: Cell index along `y`.
        :type j: int
        :returns: The Cartesian coordinate for the cell center.
        :return type: numpy.array

        See also: :func:`getSouthFaceCoordinate()`,
        :func:`getNorthFaceCoordinate()`, :func:`getEastFaceCoordinate()`,
        :func:`getWestFaceCoordinate()`, :func:`getSouthWestCornerCoordinate()`,
        :func:`getSouthEastCornerCoordinate()`,
        :func:`getNorthWestCornerCoordinate()`, and
        :func:`getNorthEastCornerCoordinate()`.
        '''
        return numpy.array([
            0.5*self.dx + i*self.dx,
            0.5*self.dy + j*self.dy]) + self.origo

    def getWestFaceCoordinate(self, i, j):
        '''
        Returns the west-oriented face-center coordinate for the cell with index
        `i` and `j`, sometimes referred to as the u-point.

        :param i: Cell index along `x`.
        :type i: int
        :param j: Cell index along `y`.
        :type j: int
        :returns: The Cartesian coordinate for the center of the western cell
            face.
        :return type: numpy.array

        See also: :func:`getCenterCoordinate()`,
        :func:`getSouthFaceCoordinate()`, :func:`getNorthFaceCoordinate()`,
        :func:`getEastFaceCoordinate()`, :func:`getSouthEastCornerCoordinate()`,
        :func:`getNorthWestCornerCoordinate()`, and
        :func:`getNorthEastCornerCoordinate()`.
        '''
        return numpy.array([
            i*self.dx,
            0.5*self.dy + j*self.dy]) + self.origo

    def getEastFaceCoordinate(self, i, j):
        '''
        Returns the east-oriented face-center coordinate for the cell with index
        `i` and `j`, sometimes referred to as the u-point.

        :param i: Cell index along `x`.
        :type i: int
        :param j: Cell index along `y`.
        :type j: int
        :returns: The Cartesian coordinate for the center of the eastern cell
            face.
        :return type: numpy.array

        See also: :func:`getCenterCoordinate()`,
        :func:`getSouthFaceCoordinate()`, :func:`getNorthFaceCoordinate()`,
        :func:`getWestFaceCoordinate()`, :func:`getSouthWestCornerCoordinate()`,
        :func:`getSouthEastCornerCoordinate()`,
        :func:`getNorthWestCornerCoordinate()`, and
        :func:`getNorthEastCornerCoordinate()`.
        '''
        return numpy.array([
            (i + 1)*self.dx,
            0.5*self.dy + j*self.dy]) + self.origo

    def getSouthFaceCoordinate(self, i, j):
        '''
        Returns the south-oriented face-center coordinate for the cell with
        index `i` and `j`, sometimes referred to as the v-point.

        :param i: Cell index along `x`.
        :type i: int
        :param j: Cell index along `y`.
        :type j: int
        :returns: The Cartesian coordinate for the center of the southern cell
            face.
        :return type: numpy.array

        See also: :func:`getCenterCoordinate()`,
        :func:`getNorthFaceCoordinate()`, :func:`getEastFaceCoordinate()`,
        :func:`getWestFaceCoordinate()`, :func:`getSouthWestCornerCoordinate()`,
        :func:`getSouthEastCornerCoordinate()`,
        :func:`getNorthWestCornerCoordinate()`, and
        :func:`getNorthEastCornerCoordinate()`.
        '''
        return numpy.array([
            0.5*self.dx + (i + 1)*self.dx,
            j*self.dy]) + self.origo

    def getNorthFaceCoordinate(self, i, j):
        '''
        Returns the north-oriented face-center coordinate for the cell with
        index `i` and `j`, sometimes referred to as the v-point.

        :param i: Cell index along `x`.
        :type i: int
        :param j: Cell index along `y`.
        :type j: int
        :returns: The Cartesian coordinate for the center of the northern cell
            face.
        :return type: numpy.array

        See also: :func:`getCenterCoordinate()`,
        :func:`getSouthFaceCoordinate()`, :func:`getEastFaceCoordinate()`,
        :func:`getWestFaceCoordinate()`, :func:`getSouthWestCornerCoordinate()`,
        :func:`getSouthEastCornerCoordinate()`,
        :func:`getNorthWestCornerCoordinate()`, and
        :func:`getNorthEastCornerCoordinate()`.
        '''
        return numpy.array([
            0.5*self.dx + (i + 1)*self.dx,
            (j + 1)*self.dy]) + self.origo

    def getSouthWestCornerCoordinate(self, i, j):
        '''
        Returns the south-west oriented corner coordinate the cell with index
        `i` and `j`, sometimes referred to as the q-point.

        :param i: Cell index along `x`.
        :type i: int
        :param j: Cell index along `y`.
        :type j: int
        :returns: The Cartesian coordinate for the center of the south-western
            cell corner.
        :return type: numpy.array

        See also: :func:`getCenterCoordinate()`,
        :func:`getSouthFaceCoordinate()`, :func:`getNorthFaceCoordinate()`,
        :func:`getEastFaceCoordinate()`, :func:`getWestFaceCoordinate()`,
        :func:`getSouthEastCornerCoordinate()`,
        :func:`getNorthWestCornerCoordinate()`, and
        :func:`getNorthEastCornerCoordinate()`.
        '''
        return numpy.array([
            i*self.dx,
            j*self.dy]) + self.origo

    def getNorthWestCornerCoordinate(self, i, j):
        '''
        Returns the north-west oriented corner coordinate the cell with index
        `i` and `j`, sometimes referred to as the q-point.

        :param i: Cell index along `x`.
        :type i: int
        :param j: Cell index along `y`.
        :type j: int
        :returns: The Cartesian coordinate for the center of the north-western
            cell corner.
        :return type: numpy.array

        See also: :func:`getCenterCoordinate()`,
        :func:`getSouthFaceCoordinate()`, :func:`getNorthFaceCoordinate()`,
        :func:`getEastFaceCoordinate()`, :func:`getWestFaceCoordinate()`,
        :func:`getSouthWestCornerCoordinate()`,
        :func:`getSouthEastCornerCoordinate()`, and
        :func:`getNorthEastCornerCoordinate()`.
        '''
        return numpy.array([
            i*self.dx,
            (j + 1)*self.dy]) + self.origo

    def getSouthEastCornerCoordinate(self, i, j):
        '''
        Returns the south-east oriented corner coordinate the cell with index
        `i` and `j`, sometimes referred to as the q-point.

        :param i: Cell index along `x`.
        :type i: int
        :param j: Cell index along `y`.
        :type j: int
        :returns: The Cartesian coordinate for the center of the south-eastern
            cell corner.
        :return type: numpy.array

        See also: :func:`getCenterCoordinate()`,
        :func:`getSouthFaceCoordinate()`, :func:`getNorthFaceCoordinate()`,
        :func:`getEastFaceCoordinate()`, :func:`getWestFaceCoordinate()`,
        :func:`getSouthWestCornerCoordinate()`,
        :func:`getNorthWestCornerCoordinate()`, and
        :func:`getNorthEastCornerCoordinate()`.
        '''
        return numpy.array([
            (i + 1)*self.dx,
            self.dy]) + self.origo

    def getNorthEastCornerCoordinate(self, i, j):
        '''
        Returns the north-east oriented corner coordinate the cell with index
        `i` and `j`, sometimes referred to as the q-point.

        :param i: Cell index along `x`.
        :type i: int
        :param j: Cell index along `y`.
        :type j: int
        :returns: The Cartesian coordinate for the center of the north-eastern
            cell corner.
        :return type: numpy.array

        See also: :func:`getCenterCoordinate()`,
        :func:`getSouthFaceCoordinate()`, :func:`getNorthFaceCoordinate()`,
        :func:`getEastFaceCoordinate()`, :func:`getWestFaceCoordinate()`,
        :func:`getSouthWestCornerCoordinate()`,
        :func:`getSouthEastCornerCoordinate()`, and
        :func:`getNorthWestCornerCoordinate()`.
        '''
        return numpy.array([
            (i + 1)*self.dx,
            (j + 1)*self.dy]) + self.origo

    def addFloe(self, icefloe):
        '''
        Add an IceFloe object to the grid.

        :param icefloe: The icefloe object to add to the grid.
        :type icefloe: seaice.IceFloeCylindrical
        '''
        if isinstance(icefloe, seaice.IceFloeCylindrical):
            self.floes.append(icefloe)
        else:
            raise Exception('Error: Incompatible icefloe data type')
