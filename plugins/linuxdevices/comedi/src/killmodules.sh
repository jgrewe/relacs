#! /bin/bash
rmmod ni_pcimio
#rmmod -f ni_mio_cs
sleep 0.1
rmmod comedi_fc
sleep 0.1
rmmod 8255
sleep 0.1
rmmod ni_tiocmd
sleep 0.1
rmmod ni_tio
sleep 0.1
rmmod mite
sleep 0.1
rmmod comedi