#! /bin/bash
rmmod dynclamp
sleep 0.1
rmmod ni_pcimio
sleep 0.1
rmmod comedi_fc
sleep 0.1
rmmod 8255
sleep 0.1
rmmod ni_tio
sleep 0.1
rmmod mite
sleep 0.1
rmmod kcomedilib
sleep 0.1
rmmod comedi
sleep 0.1
rmmod rtai_fifos.ko
sleep 0.1
rmmod rtai_up
sleep 0.1
rmmod rtai_hal
sleep 0.1
rmmod rtai_sem

