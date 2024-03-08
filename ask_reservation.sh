#!/bin/bash

# Prompt the user for a reservation
read -p "Enter your reservation: " reservation

# Perform the srun command with the reservation
srun --reservation=ST76S4 --exclusive -N 1 --pty bash