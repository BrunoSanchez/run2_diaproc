#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  copy_subtraction_files.py
#  
#  Copyright 2019 bruno <bruno.sanchez@duke.edu>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
# NEEDS TO RUN USING desc-dia 

import os
from lsst.daf.persistence import Butler


def main(visit, detector, diarepo, outputdir='.'):
    diabutler = Butler(diarepo)

    diff_img_path = os.path.join(outputdir, f'diff_exposure_v{visit}_d{detector}.fits')
    science_img_path = os.path.join(outputdir, f'science_exposure_v{visit}_d{detector}.fits')

    imgD = diabutler.get('deepDiff_differenceExp', visit=visit, detector=detector)
    imgD.writeFits(diff_img_path)

    imgS = diabutler.get('calexp', visit=visit, detector=detector)
    imgS.writeFits(science_img_path)

    print(f'files copied to {outputdir}')
    return


if __name__=='__main__':
    import argparse
    DESC = """This script will query the butler interface looking for the images
    used for DIA analysis and copy them into a specified location. This tool is 
    a simple way to compare images directly into DS9 and is not meant to be the 
    preferred way to check images.
    """
    EPIL = """Will copy images to desired location for particular visit and detector"""
    parser = argparse.ArgumentParser(description=DESC, epilog=EPIL)

    parser.add_argument('--visit', type=int, help='Visit number')
    parser.add_argument('--detector', type=int, help='Detector number')
    parser.add_argument('--diarepo', type=str, help='Data repo path')
    parser.add_argument('--outputdir', type=str, help='Output directory path', default='.')
    
    args = parser.parse_args()

    main(visit=args.visit, detector=args.detector, diarepo=args.diarepo, 
         outputdir=args.outputdir)
