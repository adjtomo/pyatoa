import os
path_pyatoa = os.path.dirname(os.path.realpath(__file__))
req_file = os.path.join(path_pyatoa, "requirements.txt")
if os.path.exists(req_file):
    with open(req_file, "r") as f:
        install_requires = list(f.read().splitlines())
else:
    install_requires = []
print(install_requires)
