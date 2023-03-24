#!/bin/bash

cp -r ../jags-tddm /tmp/

cd /tmp/jags-tddm

autoreconf -fvi && ./configure --prefix=/usr

make && sudo make install

sudo cp /usr/lib/JAGS/modules-4/tddm.* /usr/lib/x86_64-linux-gnu/JAGS/modules-4/ #todo
