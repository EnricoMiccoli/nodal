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
    def test_doc_1_6_1(self):
        path = "doc/1.6.1.csv"
        expected = (
            "Ground node: g\n"
            "e(1) \t= 2.0\n"
            "e(2) \t= -1.0\n"
            "e(4) \t= 8.0\n"
            "i(d1) \t= -1.9999999999999998\n"
            "i(e1) \t= 3.0\n"
        )
        self.assert_print(path, expected)

    def test_doc_buffer(self):
        path = "doc/buffer.csv"
        expected = (
            "Ground node: g\n"
            "e(1) \t= 9.999900000999991\n"
            "e(2) \t= 9.999900000899993\n"
            "e(3) \t= 10.0\n"
            "i(d1) \t= -9.999889805101247e-12\n"
            "i(vs) \t= 9.999900000899993e-12\n"
        )
        self.assert_print(path, expected)

    def test_doc_netlist(self):
        path = "doc/netlist.csv"
        expected = "Ground node: 1\n" "e(2) \t= -1.0\n" "e(3) \t= -2.0\n"
        self.assert_print(path, expected)

    def test_doc_opmodel_amplifier(self):
        path = "doc/opmodel_amplifier.csv"
        expected = (
            "Ground node: g\n"
            "e(1) \t= 0.9998800143982738\n"
            "e(2) \t= 1.9997600287845492\n"
            "e(3) \t= 1.0000000000000002\n"
            "e(q1_internal_node) \t= 11.998560172647306\n"
            "i(q1_vcvs) \t= 0.9998800143862756\n"
            "i(v1) \t= 1.1998560172647305e-11\n"
        )
        self.assert_print(path, expected)

    def test_doc_opmodel_voltage_buffer(self):
        path = "doc/opmodel_voltage_buffer.csv"
        expected = (
            "Ground node: g\n"
            "e(2) \t= 0.999990000099999\n"
            "e(3) \t= 0.9999999999999999\n"
            "e(q1_internal_node) \t= 0.9999900000899992\n"
            "i(q1_vcvs) \t= -9.999917560676863e-13\n"
            "i(v1) \t= 9.999900000899992e-13\n"
        )
        self.assert_print(path, expected)


class InputTesters(unittest.TestCase):
    def test_check_input_component(self):
        bad_inputs = [
            "aaaaa",  # not enough arguments
            "v1,VCVS,5,1,2",
            "v1,VCCS,5,1,2",
            "v1,CCVS,5,1,2",
            "v1,CCCS,5,1,2",
            "q1,OPMODEL,0,2,g,3",
            "v1,VCVS,5,1,2,1,1,1",  # too many arguments
            "r1,R,5,1,2,3",
            "r1,A,5,1,2,3",
            "r1,E,5,1,2,3",
            "q1,OPMODEL,1,2,g,3,1,5",
            "v1,VoltageSource,5,1,2",  # unknown type
            "r1,R,one_ohm,1,2",  # not a float
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
            "q1,OPMODEL,1,2,g,3,1",
            "q1,OPMODEL,0,2,g,3,2",
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

    def test_check_input_component_empty_line(self):
        try:
            n.check_input_component([])
        except ValueError:
            assert False

    def test_check_input_component_comment(self):
        good = "# This is a comment"
        try:
            n.check_input_component(good)
        except ValueError:
            assert False


class GroundNode(unittest.TestCase):
    def test_explicit_ground_node(self):
        inputs = [{"g": 1}, {"g": 1, "a": 10, "b": 2}]
        for deg in inputs:
            assert n.find_ground_node(deg) == "g"

    def test_highest_degree_ground_node(self):
        inputs = [({"1": 1}, "1"), ({"3": 1, "a": 10, "b": 2}, "a")]
        for deg, expected in inputs:
            assert n.find_ground_node(deg) == expected

    def test_stable_sorting_ground_node(self):
        inputs = [({"1": 1, "2": 1}, "1"), ({"3": 1, "a": 10, "b": 10}, "a")]
        for deg, expected in inputs:
            assert n.find_ground_node(deg) == expected


if __name__ == "__main__":
    logging.getLogger().setLevel(logging.CRITICAL)
    unittest.main()
