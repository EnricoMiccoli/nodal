import io
import logging
import unittest
import unittest.mock

import nodal as n


class IntegratedTest(unittest.TestCase):
    @unittest.mock.patch("sys.stdout", new_callable=io.StringIO)
    def assert_print(self, path, expected, mock_stdout):
        netlist = n.Netlist(path)
        circuit = n.Circuit(netlist)
        solution = circuit.solve()
        print(solution)
        output = mock_stdout.getvalue()
        self.assertEqual(output, expected)

    # Triple caution!!!
    # For mock.patch to work correctly, test cases must be named test_something.
    # I find this horrible, might change later
    def test_print(self):
        path = "doc/1.6.1.csv"
        expected = "Ground node: g\ne(1) \t= 2.0\ne(2) \t= -1.0\ne(4) \t= 8.0\ni(d1) \t= -1.9999999999999998\ni(e1) \t= 3.0\n"
        self.assert_print(path, expected)


class InputTesters(unittest.TestCase):
    def test_check_input_component(self):
        bad_inputs = [
            "aaaaa",  # not enough arguments
            "v1,VCVS,5,1,2",
            "v1,VCCS,5,1,2",
            "v1,CCVS,5,1,2",
            "v1,CCCS,5,1,2",
            "v1,VCVS,5,1,2,1,1,1",  # too many arguments
            "r1,R,5,1,2,3",
            "r1,A,5,1,2,3",
            "r1,E,5,1,2,3",
            "v1,VoltageSource,5,1,2",  # unknown type
            # TODO 'r1,R,one_ohm,1,2', # not a float
        ]
        good_inputs = [
            "r1,R,2,1,4",
            "r2,R,2,1,g",
            "r3,R,0.5,1,2",
            "e1,E,8,4,g",
            "a1,A,4,1,2",
            "d1,CCCS,2,2,g,1,g,r2",
            "Ri,R,1e7,1,3",
            "Ro,R,1e1,1,2",
            "vs,E,10,3,g",
            "d1,VCVS,1e5,2,g,3,1",
            "Ri,R,1e7,1,3",
            "Ro,R,1e1,1,2",
            "vs,E,10,3,g",
            "d1,VCVS,1e5,2,g,3,1",
            # TODO '' # skip empty line
        ]

        for bad in bad_inputs:
            try:
                n.check_input_component(bad.split(","))
            except ValueError:
                continue
            assert False

        for good in good_inputs:
            try:
                n.check_input_component(good.split(","))
            except ValueError:
                assert False


if __name__ == "__main__":
    logging.getLogger().setLevel(logging.CRITICAL)
    unittest.main()
