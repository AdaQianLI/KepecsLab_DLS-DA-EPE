/*******************************************************************************
 *** Filename:    MotorShield.h                                              ***
 *** Created on:  2019-11-30                                                 ***
 *** Created by:  Michael Wulf                                               ***
 ***              Cold Spring Harbor Laboratory                              ***
 ***              Kepecs Lab                                                 ***
 ***              wulf@cshl.edu                                              ***
 *** License:     GNU GPLv3                                                  ***
 ***                                                                         ***
 *** Version:     1.0.0                                                      ***
 *** Updated on:  2019-11-30                                                 ***
 *** Info:        Arduino compatible library to interface with the Arduino   ***
 ***              Motor Shield R3.                                           ***
 ***                                                                         ***
 ******************************************************************************/
/* Disclaimer:
    Arduino compatible library to interface with the Arduino Motor Shield R3.
    Copyright (C) 2019 Michael Wulf, Cold Spring Harbor Laboratory, 
    Cold Spring Harbor, New York, USA
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
  */
#ifndef MOTORSHIELD_H_
#define MOTORSHIELD_H_

/*******************************************************************************
 * Necessary includes
 ******************************************************************************/
#if defined(ARDUINO) && ARDUINO >= 100
#include "Arduino.h"
#else
#include "WProgram.h"
#endif

#include "MotorShieldDefs.h"

// Global Defines
#define MOTOR_A 0x01
#define MOTOR_B 0x02
#define MOTOR_AB 0x03

class MotorShield
{
public:
    MotorShield();
    void begin(uint8_t motor);
    void set_forward_direction(uint8_t motor, uint8_t direction);
    void brake(uint8_t motor);
    void stop(uint8_t motor);
    void forward(uint8_t motor, uint8_t speed);
    void backward(uint8_t motor, uint8_t speed);

private:
    uint8_t _forward_dir_A;
    uint8_t _forward_dir_B;
};
#endif