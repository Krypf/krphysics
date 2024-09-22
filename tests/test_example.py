import unittest

def hello_world():
    return "Hello, World!"

class TestExampleModule(unittest.TestCase):
    def test_hello_world(self):
        self.assertEqual(hello_world(), "Hello, World!")

if __name__ == "__main__":
    unittest.main()
