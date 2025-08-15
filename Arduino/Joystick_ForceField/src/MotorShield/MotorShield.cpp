/*******************************************************************************
 *** Filename:    MotorShield.cpp                                            ***
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

#include "MotorShield.h"

// Constructor
MotorShield::MotorShield()
{
    this->_forward_dir_A = HIGH;
    this->_forward_dir_B = HIGH;
}

// Initialize the object
void MotorShield::begin(uint8_t motor)
{
    if (motor & MOTOR_A)
    {
        // GPIOs for motor A
        pinMode(MOTOR_A_DIR_PIN, OUTPUT);
        pinMode(MOTOR_A_BRK_PIN, OUTPUT);
        pinMode(MOTOR_A_PWM_PIN, OUTPUT);

        // Set PWM output to 0 to stop motors
        analogWrite(MOTOR_A_PWM_PIN, 0);
    }
    if (motor & MOTOR_B)
    {
        // GPIOs for motor B
        pinMode(MOTOR_B_DIR_PIN, OUTPUT);
        pinMode(MOTOR_B_BRK_PIN, OUTPUT);
        pinMode(MOTOR_B_PWM_PIN, OUTPUT);

        // Set PWM output to 0 to stop motors
        analogWrite(MOTOR_B_PWM_PIN, 0);
    }
}

// Change the direction for a given motor
void MotorShield::set_forward_direction(uint8_t motor, uint8_t direction)
{
    if (motor & MOTOR_A)
    {
        this->_forward_dir_A = !(!direction);
    }
    if (motor & MOTOR_B)
    {
        this->_forward_dir_B = !(!direction);
    }
}

void MotorShield::brake(uint8_t motor)
{
    if (motor & MOTOR_A)
    {
        // Set PWM output to 0
        analogWrite(MOTOR_A_PWM_PIN, 0);
        // Set Brake pin to HIGH
        digitalWrite(MOTOR_A_BRK_PIN, HIGH);
    }
    if (motor & MOTOR_B)
    {
        // Set PWM output to 0
        analogWrite(MOTOR_B_PWM_PIN, 0);
        // Set Brake pin to HIGH
        digitalWrite(MOTOR_B_BRK_PIN, HIGH);
    }
}

void MotorShield::stop(uint8_t motor)
{
    if (motor & MOTOR_A)
    {
        // Set PWM output to 0
        analogWrite(MOTOR_A_PWM_PIN, 0);
    }
    if (motor & MOTOR_B)
    {
        // Set PWM output to 0
        analogWrite(MOTOR_B_PWM_PIN, 0);
    }
}

void MotorShield::forward(uint8_t motor, uint8_t speed)
{
    if (motor & MOTOR_A)
    {
        // Disable brake
        digitalWrite(MOTOR_A_BRK_PIN, LOW);

        // Set direction
        digitalWrite(MOTOR_A_DIR_PIN, (this->_forward_dir_A));

        // Set PWM output
        analogWrite(MOTOR_A_PWM_PIN, speed);
    }
    if (motor & MOTOR_B)
    {
        // Disable brake
        digitalWrite(MOTOR_B_BRK_PIN, LOW);

        // Set direction
        digitalWrite(MOTOR_B_DIR_PIN, (this->_forward_dir_B));

        // Set PWM output
        analogWrite(MOTOR_B_PWM_PIN, speed);
    }
}

void MotorShield::backward(uint8_t motor, uint8_t speed)
{
    if (motor & MOTOR_A)
    {
        // Disable brake
        digitalWrite(MOTOR_A_BRK_PIN, LOW);

        // Set direction
        digitalWrite(MOTOR_A_DIR_PIN, !(this->_forward_dir_A));

        // Set PWM output
        analogWrite(MOTOR_A_PWM_PIN, speed);
    }
    if (motor & MOTOR_B)
    {
        // Disable brake
        digitalWrite(MOTOR_B_BRK_PIN, LOW);

        // Set direction
        digitalWrite(MOTOR_B_DIR_PIN, !(this->_forward_dir_B));

        // Set PWM output
        analogWrite(MOTOR_B_PWM_PIN, speed);
    }
}