import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import matplotlib.patches as pat
from scipy import stats
import platform
import os
import matplotlib

if __name__ == '__main__':
    #os.chdir('../')
    path = os.getcwd()
    opt_data = pd.read_csv(path+"/build/res.csv")
    obs_data = pd.read_csv(path+"/build/obs.csv")

    print(opt_data)
    print(obs_data.columns)
    

    linewidth_ = 5.5

    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['mathtext.fontset'] = 'stix' # math fontの設定
    plt.rcParams['xtick.direction'] = 'in' # x axis in
    plt.rcParams['ytick.direction'] = 'in' # y axis in
    plt.rcParams['axes.linewidth'] = 1.0 # axis line width
    plt.rcParams['axes.grid'] = True # make grid
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    fig = plt.figure(figsize=(80, 40))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    #ax1.plot(obs_data['obs_time'], obs_data['obs_s'], label="Obstacle", color="black", linewidth=linewidth_, linestyle="dashed")
    ax1.plot(opt_data['time'], opt_data['original_dis'], label="Original Dis", color="black", linewidth=linewidth_)
    ax1.plot(opt_data['time'], opt_data['op_dis'], label="Optimized Dis", color="purple", linewidth=linewidth_)
    #ax1.plot(opt_data['jerk_filtered_time'], opt_data['position'], label="Jerk Filter", color="orange", linewidth=linewidth_)
    #ax1.plot(opt_data['lp_time'], opt_data['lp_position'], label="LP", color="blue", linewidth=linewidth_)
    if (abs(obs_data['s_front'][0] - obs_data['s_front'][1])> 0.01): # exist
        ax1.plot(obs_data['time'], obs_data['s_front'], label="front_veh", color="red", linewidth=0.5*linewidth_)
        ax1.plot(obs_data['time'], obs_data['s_front']-4, color="red", linewidth=0.5*linewidth_)
    if (abs(obs_data['s_behind'][0] - obs_data['s_behind'][1])> 0.01): # exist
        ax1.plot(obs_data['time'], obs_data['s_behind'], label="behind_veh", color="green", linewidth=0.5*linewidth_)
        ax1.plot(obs_data['time'], obs_data['s_behind']+4, color="green", linewidth=0.5*linewidth_)
    if (abs(obs_data['s_egolanefront'][0] - obs_data['s_egolanefront'][1])> 0.01): # exist
        ax1.plot(obs_data['time'], obs_data['s_egolanefront'], label="egolanefront_veh", color="blue", linewidth=0.5*linewidth_)
        ax1.plot(obs_data['time'], obs_data['s_egolanefront']-4, color="blue", linewidth=0.5*linewidth_)
    ax1.set_xlabel("t [s]", fontsize=50)
    ax1.set_ylabel("s [m]", fontsize=50)
    ax1.tick_params(labelsize=40)
    ax1.legend(fontsize=40)

    

    #ax1.plot(obs_data['obs_time'], obs_data['obs_s'], label="Obstacle", color="black", linewidth=linewidth_, linestyle="dashed")
    ax2.plot(opt_data['time'], opt_data['original_vel'], label="Original Velocity", color="black", linewidth=linewidth_)
    ax2.plot(opt_data['time'], opt_data['op_vel'], label="Optimized Velocity", color="blue", linewidth=linewidth_)
    #ax1.plot(opt_data['jerk_filtered_time'], opt_data['position'], label="Jerk Filter", color="orange", linewidth=linewidth_)
    #ax1.plot(opt_data['lp_time'], opt_data['lp_position'], label="LP", color="blue", linewidth=linewidth_)
    ax2.set_xlabel("t [s]", fontsize=50)
    ax2.set_ylabel("v [m/s]", fontsize=50)
    ax2.tick_params(labelsize=40)
    ax2.legend(fontsize=40)

    ax3.plot(opt_data['time'], opt_data['original_acc'], label="Original Acceleration", color="black", linewidth=linewidth_)
    ax3.plot(opt_data['time'], opt_data['op_acc'], label="Optimized Acceleration", color="orange", linewidth=linewidth_)
    #ax1.plot(opt_data['jerk_filtered_time'], opt_data['position'], label="Jerk Filter", color="orange", linewidth=linewidth_)
    #ax1.plot(opt_data['lp_time'], opt_data['lp_position'], label="LP", color="blue", linewidth=linewidth_)
    ax3.set_xlabel("t [s]", fontsize=50)
    ax3.set_ylabel("a [m/s2]", fontsize=50)
    ax3.tick_params(labelsize=40)
    ax3.legend(fontsize=40)

    ax4.plot(opt_data['time'], opt_data['original_jerk'], label="Original Jerk", color="black", linewidth=linewidth_)
    ax4.plot(opt_data['time'], opt_data['op_jerk'], label="Optimized Jerk", color="red", linewidth=linewidth_)
    #ax1.plot(opt_data['jerk_filtered_time'], opt_data['position'], label="Jerk Filter", color="orange", linewidth=linewidth_)
    #ax1.plot(opt_data['lp_time'], opt_data['lp_position'], label="LP", color="blue", linewidth=linewidth_)
    ax4.set_xlabel("t [s]", fontsize=50)
    ax4.set_ylabel("a [m/s3]", fontsize=50)
    ax4.tick_params(labelsize=40)
    ax4.legend(fontsize=40)


    #userhome = os.path.expanduser('~')
    #desktop_path = userhome + '/Desktop/'
    plt.savefig('experiment_obstacle_avoidance.png', bbox_inches="tight", pad_inches=0.05)