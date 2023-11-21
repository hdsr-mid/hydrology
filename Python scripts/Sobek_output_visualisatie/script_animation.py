# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 11:43:49 2023

@author: PetraH
"""

import cv2
import os
import config

if __name__ == "__main__":
    suffices = ['2u_GLG','24u_GLG','96u_GLG','2u_GHG','24u_GHG','96u_GHG']
    for suffix in suffices:    
        image_folder = os.path.join(config.output,'figs')
        video_name   = os.path.join(config.output,'video_'+suffix+'.mp4')
        
        images = ['fig_'+suffix + '_t' + str(t)+'.png' for t in range(0,151)]
        frame = cv2.imread(os.path.join(image_folder, images[0]))
        height, width, layers = frame.shape
        
        video = cv2.VideoWriter(video_name, 0, 1, (width,height))
        
        for image in images:
            video.write(cv2.imread(os.path.join(image_folder, image)))
        
        cv2.destroyAllWindows()
        video.release()