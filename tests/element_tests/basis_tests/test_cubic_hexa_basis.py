import unittest

from tacs import TACS, elements


class BasisTest(unittest.TestCase):
    def setUp(self):
        # fd/cs step size
        if TACS.dtype is complex:
            self.dh = 1e-50
            self.rtol = 1e-11
        else:
            self.dh = 1e-6
            self.rtol = 1e-2

        # Basically, only check relative tolerance
        self.atol = 1e99
        self.print_level = 0

        # Create the basis functions for 3D
        self.basis = elements.CubicHexaBasis()

        # Seed random number generator in tacs for consistent test results
        elements.SeedRandomGenerator(0)

    def test_element_basis_functions(self):
        fail = elements.TestElementBasisFunctions(
            self.basis, self.dh, self.print_level, self.atol, self.rtol
        )
        self.assertFalse(fail)

    def test_element_basis_face_normals(self):
        fail = elements.TestElementBasisFaceNormals(
            self.basis, self.dh, self.print_level, self.atol, self.rtol
        )
        self.assertFalse(fail)

    def test_element_basis_jacobian_transform(self):
        fail = elements.TestElementBasisJacobianTransform(
            self.basis, self.dh, self.print_level, self.atol, self.rtol
        )
        self.assertFalse(fail)
