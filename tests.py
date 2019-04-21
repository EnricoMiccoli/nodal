import unittest
import unittest.mock
import nodal as n
import io
import sys

class IntegratedTest(unittest.TestCase):
    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
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

if __name__=="__main__":
    unittest.main()

