''' This file is part of cagl.
 *
 * cagl is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * cagl is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with cagl.  If not, see <https://www.gnu.org/licenses/>
 *
 * Authors:
 * Lorenzo Magherini (m4gh3) '''

import torch
import torch.nn as nn
import torch.nn.functional as F
import math


class PGate(nn.Module):

    def __init__(self, in_c, out_c, padding=0 ):

        super().__init__()

        self.conv0  = nn.Conv2d(in_c , out_c, 1, bias=False )
        self.conv1  = nn.Conv2d(out_c, out_c, 5, stride=2, padding=padding, bias=False, groups=out_c )
        self.pxus   = nn.PixelUnshuffle(2)
        self.conv2  = nn.Conv2d(4*out_c, out_c, 1, bias=False )
        self.conv3  = nn.Conv2d(4*out_c, out_c, 1, bias=False )
        self.conv4  = nn.Conv2d(4*out_c, out_c, 1, bias=False )


        self.bn0    = nn.BatchNorm2d(out_c)
        self.bn1    = nn.BatchNorm2d(out_c)
        self.bn2    = nn.BatchNorm2d(out_c)

        
        torch.nn.init.dirac_(self.conv2.weight)
        torch.nn.init.dirac_(self.conv3.weight) 
        torch.nn.init.dirac_(self.conv4.weight)  


    def forward(self, x ):

        x  = self.pxus(self.conv1(self.conv0(x)))
        x0 = self.bn0(self.conv2(x))/2+0.5
        x1 = self.bn1(self.conv3(x))/2+0.5
        x2 = self.bn2(self.conv4(x))/2


        x = 0.1515*x0*x1 + 0.7376*x2

        return x


class NetModule(nn.Module):

    def __init__(self, conf ):
        super().__init__()

        #self.pgates = nn.ModuleList([PGate(3, 4 ), PGate(4, 4 ), PGate(4, 8 ), PGate(8, 8 )])
        self.pgates = nn.ModuleList([PGate(3, 4, 8 ), PGate(4, 2, 0)])#, PGate(8, 16 ), PGate(16, 16 ), PGate(16, 32 )])
        self.dense  = nn.Linear(8, 10 )

    def forward(self, x ):
 
        for pgate in self.pgates:
            x = pgate(x)

        x = torch.flatten(x, 1 )
        x = self.dense(x)

        return x

