#!/bin/sh

## generate locked environment file with full specification of versions

conda env export --file environment.lock.yml -n HBVouroboros
