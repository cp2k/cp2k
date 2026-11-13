#!/usr/bin/env python3

# author: Ole Schuett

import sys
import yaml


# ======================================================================================
def main() -> None:
    with open(sys.argv[1]) as f:
        doc = yaml.safe_load(f)

    for d in doc["Diagnostics"]:
        m = d["DiagnosticMessage"]
        filename = m["FilePath"].split("/cp2k/src/", 1)[1]
        message = m["Message"]
        name = d["DiagnosticName"]
        print(f"{filename}: {message} https://cp2k.org/conv#{name}")


# ======================================================================================
if __name__ == "__main__":
    main()

# EOF
