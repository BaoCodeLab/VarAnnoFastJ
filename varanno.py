#!/usr/bin/env python
### Command-line interface for VarAnnoFastJ

from annolib import parameters

if __name__ == "__main__":
    print("here")
    args=parameters.parse_arguments()
    print(args)
    args.func(args)




