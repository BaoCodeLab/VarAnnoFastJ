#!/usr/bin/env python
### Command-line interface for VarAnnoFastJ

from annolib import parameters2

if __name__ == "__main__":
    args=parameters2.parse_arguments()
    args.func(args)




