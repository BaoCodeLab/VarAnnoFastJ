#!/usr/bin/env python
### Command-line interface for VarAnnoFastJ

from annolib import parameters3

if __name__ == "__main__":
    args=parameters3.parse_arguments()
    args.func(args)




