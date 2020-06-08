#!/usr/bin/env python
### Command-line interface for VarAnnoFastJ

from annolib import parameters

if __name__ == "__main__":
    args=parameters.parse_arguments()
    args.func(args)




