#!/usr/bin/env python
import numpy


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
        self.nx = nx
        self.ny = ny

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
            self.dx = dx
            self.Lx = dx * self.nx
        else:
            self.Lx = Lx
            self.dx = Lx/float(self.nx)
        if dy:
            self.dy = dy
            self.Ly = dy * self.ny
        else:
            self.Ly = Ly
            self.dy = Ly/float(self.ny)

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
        '''
        return numpy.array([
            0.5*self.dx + i*self.dx,
            0.5*self.dy + j*self.dy])

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
        '''
        return numpy.array([
            i*self.dx,
            0.5*self.dy + j*self.dy])

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
        '''
        return numpy.array([
            (i + 1)*self.dx,
            0.5*self.dy + j*self.dy])

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
        '''
        return numpy.array([
            0.5*self.dx + (i + 1)*self.dx,
            j*self.dy])

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
        :func:`getSouthFaceCoordinate()`, :func:`getNorthFaceCoordinate()`,
        :func:`getWestFaceCoordinate()`, :func:`getWestFaceCoordinate()`,
        :func:`getSouthWestCornerCoordinate()`,
        :func:`getSouthEastCornerCoordinate()`,
        :func:`getNorthWestCornerCoordinate()`,
        :func:`getNorthEastCornerCoordinate()`
        '''
        return numpy.array([
            0.5*self.dx + (i + 1)*self.dx,
            (j + 1)*self.dy])

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
        :func:`getWestFaceCoordinate()`, :func:`getWestFaceCoordinate()`,
        :func:`getSouthWestCornerCoordinate()`,
        :func:`getSouthEastCornerCoordinate()`,
        :func:`getNorthWestCornerCoordinate()`,
        :func:`getNorthEastCornerCoordinate()`
        '''
        return numpy.array([
            i*self.dx,
            j*self.dy])

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
        :func:`getWestFaceCoordinate()`, :func:`getWestFaceCoordinate()`,
        :func:`getSouthWestCornerCoordinate()`,
        :func:`getSouthEastCornerCoordinate()`,
        :func:`getNorthWestCornerCoordinate()`,
        :func:`getNorthEastCornerCoordinate()`
        '''
        return numpy.array([
            i*self.dx,
            (j + 1)*self.dy])

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
        :func:`getWestFaceCoordinate()`, :func:`getWestFaceCoordinate()`,
        :func:`getSouthWestCornerCoordinate()`,
        :func:`getSouthEastCornerCoordinate()`,
        :func:`getNorthWestCornerCoordinate()`,
        :func:`getNorthEastCornerCoordinate()`
        '''
        return numpy.array([
            (i + 1)*self.dx,
            self.dy])

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
        :func:`getWestFaceCoordinate()`, :func:`getWestFaceCoordinate()`,
        :func:`getSouthWestCornerCoordinate()`,
        :func:`getSouthEastCornerCoordinate()`,
        :func:`getNorthWestCornerCoordinate()`,
        :func:`getNorthEastCornerCoordinate()`
        '''
        return numpy.array([
            (i + 1)*self.dx,
            (j + 1)*self.dy])
