#! /bin/sh
#
# build.sh
# Copyright (C) 2019 Malcolm Ramsay <malramsay64@gmail.com>
#
# Distributed under terms of the MIT license.
#


#!/bin/bash
cargo build --release -p crystal-packing-cli
cargo install --bin packing --path crystal-packing-cli --root $PREFIX
