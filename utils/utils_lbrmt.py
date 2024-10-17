import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib.collections import LineCollection
import os


def load_file(outdir, fr):
    ''' Load one .txt file. '''
    frname = "/fr_%04d" % fr
    filename = outdir+frname+".txt"

    contents = np.loadtxt(filename)
    header = contents[0, :]
    nx = int(header[-2])
    ny = int(header[-1])
    rho = contents[1:nx*ny+1, 2].reshape(ny, nx)
    ux = contents[1:nx*ny+1, 3].reshape(ny, nx)
    uy = contents[1:nx*ny+1, 4].reshape(ny, nx)
    speed = np.sqrt(ux*ux+uy*uy)
    vort = np.gradient(uy, axis=1)-np.gradient(ux, axis=0)
    Fx = contents[1:nx*ny+1, 5].reshape(ny, nx)
    Fy = contents[1:nx*ny+1, 6].reshape(ny, nx)
    ss = contents[1:nx*ny+1, 7].reshape(ny, nx)

    return nx, ny, rho, ux, uy, speed, vort, Fx, Fy, ss


def load_mmap_file(outdir, fr):
    ''' Load one multimaps .txt file.'''
    frname = "/fr_%04d" % fr
    filename = outdir+"/mmap"+frname+".txt"
    contents = np.loadtxt(filename)

    header = contents[0, :]
    nx = int(header[-3])
    ny = int(header[-2])
    nobjs = int(header[-1])

    # Initialize fields for individual object
    phi_dict = {}
    remapX_dict = {}
    remapY_dict = {}
    cc_dict = {}
    detF_dict = {}
    for obj_id in range(nobjs):
        phi_dict[obj_id] = np.ones((ny, nx))*nx*ny
        remapX_dict[obj_id] = np.empty((ny, nx))
        remapY_dict[obj_id] = np.empty((ny, nx))
        cc_dict[obj_id] = -np.ones((ny, nx))
        detF_dict[obj_id] = np.empty((ny, nx))

    for row in range(1, contents.shape[0]):
        obj_i = int(contents[row][1])
        obj_j = int(contents[row][2])
        obj_remapX = contents[row][3]
        obj_remapY = contents[row][4]
        obj_phi = contents[row][5]
        obj_cc = contents[row][6]
        obj_detF = contents[row][7]

        # Each object as dictionary value
        obj_id = int(contents[row][0])
        for nobj in range(nobjs):
            if obj_id == nobj:
                remapX_dict[obj_id][obj_j, obj_i] = obj_remapX
                remapY_dict[obj_id][obj_j, obj_i] = obj_remapY
                phi_dict[obj_id][obj_j, obj_i] = obj_phi
                cc_dict[obj_id][obj_j, obj_i] = obj_cc
                detF_dict[obj_id][obj_j, obj_i] = obj_detF

    return nobjs, remapX_dict, remapY_dict, phi_dict, cc_dict, detF_dict


def load_mmap_file_full(outdir, fr):
    ''' Load one multimaps .txt file in one np array.'''
    filename = "/fr_%04d" % fr
    contents = np.loadtxt(outdir+"/mmap"+filename+".txt")
    header = contents[0, :]
    nx = int(header[-3])
    ny = int(header[-2])

    # Initialize fields for all objects
    phi = np.ones((ny, nx))*nx*ny
    remapX = np.empty((ny, nx))
    remapY = np.empty((ny, nx))
    cc = -np.ones((ny, nx))

    for row in range(1, contents.shape[0]):
        obj_i = int(contents[row][1])
        obj_j = int(contents[row][2])
        obj_remapX = contents[row][3]
        obj_remapY = contents[row][4]
        obj_phi = contents[row][5]
        obj_cc = contents[row][6]
        # Ensure solid remap won't be overwritten by
        # extrapolation remap of other solids
        if obj_phi < 1.:
            remapX[obj_j, obj_i] = obj_remapX
            remapY[obj_j, obj_i] = obj_remapY
            phi[obj_j, obj_i] = obj_phi
            cc[obj_j, obj_i] = obj_cc

    return remapX, remapY, phi, cc


def create_dirs(outdir, save_plots=False, save_results=False):
    ''' Create directories for plots. '''
    if save_plots:
        plots_dir = os.path.join(outdir, 'plots')
    if save_results:
        dirname = os.path.dirname(outdir)
        plots_dir = os.path.join(dirname, 'results')
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
    return plots_dir


def colorbar_limits(outdir, fr_end):
    ''' Create limits for color bars based on the last frame. '''
    # Density
    rho = load_file(outdir, fr_end)[2]
    rho_diff = max(1.-np.min(rho), np.max(rho)-1.)
    rho_min = 1.-rho_diff
    rho_max = 1.+rho_diff
    # rho_min=np.min(rho)
    # rho_max=np.max(rho)
    # Velocity
    speed = load_file(outdir, fr_end)[5]
    speed_min = np.min(speed)
    speed_max = np.max(speed)
    # Vorticity
    vort = load_file(outdir, fr_end)[6]
    vort_min = -max(abs(np.min(vort)), abs(np.max(vort)))
    vort_max = -vort_min

    return rho_min, rho_max, speed_min, speed_max, vort_min, vort_max


def max_limits(field, outdir, fr_end):
    ''' Create limits for the color bars based on the max value of the field. '''
    if field == 'stress':
        ss_max = []
        for fr in range(fr_end+1):
            ss = load_file(outdir, fr)[9]
            ss_max.append(np.max(ss))
        return max(ss_max)


def plot_field(field, outdir, fr, dt, figsize, dpi, cmap, cmin, cmax, cvals,
               plot_rmap=False, plot_streamlines=False,
               plot_cbar=False, plot_title=False, plot_ticks=True,
               save_plots=False, save_results=False):
    if field == 'density':
        return plot_density(outdir, fr, dt, figsize, dpi, cmap, cmin, cmax, cvals,
                            plot_rmap, plot_streamlines, plot_cbar, plot_title, plot_ticks, save_plots, save_results)
    elif field == 'velocity':
        return plot_velocity(outdir, fr, dt, figsize, dpi, cmap, cmin, cmax, cvals,
                             plot_rmap, plot_streamlines, plot_cbar, plot_title, plot_ticks, save_plots, save_results)
    elif field == 'vorticity':
        return plot_vorticity(outdir, fr, dt, figsize, dpi, cmap, cmin, cmax, cvals,
                              plot_rmap, plot_streamlines, plot_cbar, plot_title, plot_ticks, save_plots, save_results)
    elif field == 'stress':
        return plot_stress(outdir, fr, dt, figsize, dpi, cmap, cmin, cmax, cvals,
                           plot_rmap, plot_streamlines, plot_cbar, plot_title, plot_ticks, save_plots, save_results)
    else:
        raise NotImplementedError(
            "Only can plot density/velocity/vorticity field!")


def plot_density(outdir, fr, dt, figsize, dpi, cmap, cmin, cmax, cvals,
                 plot_rmap=False, plot_streamlines=False,
                 plot_cbar=False, plot_title=False, plot_ticks=True,
                 save_plots=False, save_results=False):
    # Load in the physical fields at frame fr
    nx, ny, rho, ux, uy, speed, vort, Fx, Fy, ss = load_file(outdir, fr)
    # Create mesh
    x = np.arange(0.5, nx+0.5, 1)
    y = np.arange(0.5, ny+0.5, 1)
    X, Y = np.meshgrid(x, y)

    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)
    # Fix axis location
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0,
                    chartBox.width,
                    chartBox.height])
    # Plot density
    img = ax.imshow(rho, cmap=cmap, interpolation='none', vmin=cmin, vmax=cmax,
                    origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])
    # Colorbar
    if plot_cbar:
        cbar = fig.colorbar(img, ax=ax, format="%.4f", pad=0.02, ticks=cvals)
        cbar.ax.set_yticklabels(cvals)

    # Load in the solid at frame fr
    nobjs, remapX_dict, remapY_dict, phi_dict, cc_dict, detF_dict = load_mmap_file(
        outdir, fr)
    # Masked velocity
    mux, muy = ux, uy
    # Reference map contour interval
    xlevels = np.arange(0, nx, 10)
    ylevels = np.arange(0, ny, 10)
    # Plot individual solid
    for n in range(nobjs):
        # Plot solid–fluid interface
        ax.contour(X, Y, phi_dict[n], levels=[
                   0.], colors='k', linewidths=0.5, linestyles='solid')
        if plot_rmap:
            # Create reference map mask
            remap_mask = (phi_dict[n] > 0.).astype(int)
            mask = (remap_mask).astype(bool)
            remapX = np.ma.array(remapX_dict[n], mask=mask)
            remapY = np.ma.array(remapY_dict[n], mask=mask)
            # Create velocity mask
            speed_mask = (phi_dict[n] < 0.).astype(int)
            mask = (speed_mask).astype(bool)
            mux = np.ma.array(mux, mask=mask)
            muy = np.ma.array(muy, mask=mask)
            # Plot reference map contours
            ax.contour(X, Y, remapX, levels=xlevels, colors='k',
                       linewidths=0.25, linestyles='solid', corner_mask=True)
            ax.contour(X, Y, remapY, levels=ylevels, colors='k',
                       linewidths=0.25, linestyles='solid', corner_mask=True)

    if plot_streamlines:
        # Overlay streamlines
        lw = speed/speed.max()*0.5
        ax.streamplot(X, Y, mux, muy, density=1.0, linewidth=lw,
                      color='#3399ff', arrowsize=0.1, arrowstyle='->', broken_streamlines=False, zorder=0)

    # Layout bookkeeping
    if plot_title:
        ax.set_title(('density'+' at '+r'$t=%i \Delta t$') % (fr*dt))

    # Set axis limits
    ax.set_aspect('equal')
    ax.set_xlim(0, nx)
    ax.set_ylim(0, ny)
    # Set axis ticks
    if plot_ticks:
        ax.set_xticks([0, nx//2, nx])
        ax.set_xticklabels([0, r'$L/2$', '$L$'], fontsize=8)
        ax.set_yticks([0, ny//2, ny])
        ax.set_yticklabels([0, r'$\dfrac{H}{2}$', '$H$'], fontsize=8)
    else:
        ax.set_xticks(np.arange(0, nx, 50))
        ax.set_xticklabels([])
        ax.set_yticks(np.arange(0, ny, 50))
        ax.set_yticklabels([])
    ax.minorticks_off()
    ax.tick_params(axis='x', which='both', bottom=True,
                   top=False, labelbottom=True, size=2)
    ax.tick_params(axis='y', which='both', left=True,
                   right=False, labelleft=True, size=2)
    ax.margins(0.5, 0.5)

    if (save_plots or save_results) and outdir != None:
        savedir = create_dirs(outdir, save_plots, save_results)
        if save_plots:
            figname = "/rho_frame_%i.png" % (fr)
            fig.savefig(savedir+figname, transparent=True, bbox_inches='tight',
                        pad_inches=0)
        if save_results:
            figname = "/" + \
                os.path.split(outdir)[-1][:-4] + "_rho_frame_%i.png" % (fr)
            fig.savefig(savedir+figname, transparent=True, bbox_inches='tight',
                        pad_inches=0.005)
        plt.close()

    return fig, ax, img


def plot_velocity(outdir, fr, dt, figsize, dpi, cmap, cmin, cmax, cvals,
                  plot_rmap=False, plot_streamlines=False,
                  plot_cbar=False, plot_title=False, plot_ticks=True,
                  save_plots=False, save_results=False):
    # Load in the physical fields at frame fr
    nx, ny, rho, ux, uy, speed, vort, Fx, Fy, ss = load_file(outdir, fr)
    # Create mesh
    x = np.arange(0.5, nx+0.5, 1)
    y = np.arange(0.5, ny+0.5, 1)
    X, Y = np.meshgrid(x, y)

    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)
    # Fix axis location
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0,
                    chartBox.width,
                    chartBox.height])
    # Place holder white background
    img = ax.imshow(np.ones(speed.shape), cmap=cm.binary,
                    origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])

    # Load in the solid at frame fr
    nobjs, remapX_dict, remapY_dict, phi_dict, cc_dict, detF_dict = load_mmap_file(
        outdir, fr)
    # Masked velocity
    mux, muy = ux, uy
    # Reference map contour interval
    xlevels = np.arange(0, nx, 5)
    ylevels = np.arange(0, ny, 5)
    # Plot individual solid
    for n in range(nobjs):
        # Plot solid–fluid interface
        ax.contour(X, Y, phi_dict[n], levels=[
                   0.], colors='#93a1ad', linewidths=0.5, linestyles='solid')
        if plot_rmap:
            # Create reference map mask
            remap_mask = (phi_dict[n] > 0.).astype(int)
            mask = (remap_mask).astype(bool)
            remapX = np.ma.array(remapX_dict[n], mask=mask)
            remapY = np.ma.array(remapY_dict[n], mask=mask)
            # Create velocity mask
            speed_mask = (phi_dict[n] < 0.).astype(int)
            mask = (speed_mask).astype(bool)
            mux = np.ma.array(mux, mask=mask)
            muy = np.ma.array(muy, mask=mask)
            # Plot reference map contours
            ax.contour(X, Y, remapX, levels=xlevels, colors='k',
                       linewidths=0.25, linestyles='solid', corner_mask=True)
            ax.contour(X, Y, remapY, levels=ylevels, colors='k',
                       linewidths=0.25, linestyles='solid', corner_mask=True)


    if plot_streamlines:
        # Overlay streamlines
        lw = speed / 0.01
        ax.streamplot(X, Y, ux, uy, density=0.36, linewidth=lw,
                      color='#3399ff', arrowsize=1.0, arrowstyle='->', broken_streamlines=False, zorder=0)

    # Layout bookkeeping
    if plot_title:
        ax.set_title(('streamlines'+' at '+r'$t=%i \Delta t$') % (fr*dt))

    # Set axis limits
    ax.set_aspect('equal')
    ax.set_xlim(0, nx)
    ax.set_ylim(0, ny)
    # Set axis ticks
    if plot_ticks:
        ax.set_xticks([0, nx//2, nx])
        ax.set_xticklabels([0, r'$L/2$', '$L$'], fontsize=8)
        ax.set_yticks([0, ny//2, ny])
        ax.set_yticklabels([0, r'$\dfrac{H}{2}$', '$H$'], fontsize=8)
    else:
        ax.set_xticks([0, nx//2, nx])
        ax.set_xticklabels([])
        ax.set_yticks([0, ny//2, ny])
        ax.set_yticklabels([])
    ax.tick_params(axis='x', which='both', bottom=True,
                   top=False, labelbottom=True)
    ax.tick_params(axis='y', which='both', left=True,
                   right=False, labelleft=True)
    ax.margins(0.5, 0.5)

    if (save_plots or save_results) and outdir != None:
        savedir = create_dirs(outdir, save_plots, save_results)
        if save_plots:
            figname = "/vel_frame_%i.png" % (fr)
            fig.savefig(savedir+figname, transparent=True, bbox_inches='tight',
                        pad_inches=0)
        if save_results:
            figname = "/" + \
                os.path.split(outdir)[-1][:-4] + "_vel_frame_%i.png" % (fr)
            fig.savefig(savedir+figname, transparent=True, bbox_inches='tight',
                        pad_inches=0.005)
        plt.close()

    return fig, ax, img


def plot_vorticity(outdir, fr, dt, figsize, dpi, cmap, cmin, cmax, cvals,
                   plot_rmap=False, plot_streamlines=False,
                   plot_cbar=False, plot_title=False, plot_ticks=True,
                   save_plots=False, save_results=False):
    # Load in the physical fields at frame fr
    nx, ny, rho, ux, uy, speed, vort, Fx, Fy, ss = load_file(outdir, fr)
    # Create mesh
    x = np.arange(0.5, nx+0.5, 1)
    y = np.arange(0.5, ny+0.5, 1)
    X, Y = np.meshgrid(x, y)

    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)
    # Fix axis location
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0,
                    chartBox.width,
                    chartBox.height])
    # Plot density
    img = ax.imshow(vort, cmap=cmap, interpolation='none', vmin=cmin, vmax=cmax,
                    origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])
    # Colorbar
    if plot_cbar:
        cbar = fig.colorbar(img, ax=ax, format="%.4f", pad=0.02, ticks=cvals)
        cbar.ax.set_yticklabels(cvals)

    # Load in the solid at frame fr
    nobjs, remapX_dict, remapY_dict, phi_dict, cc_dict, detF_dict = load_mmap_file(
        outdir, fr)
    # Masked velocity
    mux, muy = ux, uy
    # Reference map contour interval
    xlevels = np.arange(0, nx, 10)
    ylevels = np.arange(0, ny, 10)
    # Plot individual solid
    for n in range(nobjs):
        # Plot solid–fluid interface
        ax.contour(X, Y, phi_dict[n], levels=[
                   0.], colors='k', linewidths=0.5, linestyles='solid')
        if plot_rmap:
            # Create reference map mask
            remap_mask = (phi_dict[n] > 0.).astype(int)
            mask = (remap_mask).astype(bool)
            remapX = np.ma.array(remapX_dict[n], mask=mask)
            remapY = np.ma.array(remapY_dict[n], mask=mask)
            # Create velocity mask
            speed_mask = (phi_dict[n] < 0.).astype(int)
            mask = (speed_mask).astype(bool)
            mux = np.ma.array(mux, mask=mask)
            muy = np.ma.array(muy, mask=mask)
            # Plot reference map contours
            ax.contour(X, Y, remapX, levels=xlevels, colors='k',
                       linewidths=0.25, linestyles='solid', corner_mask=True)
            ax.contour(X, Y, remapY, levels=ylevels, colors='k',
                       linewidths=0.25, linestyles='solid', corner_mask=True)

    if plot_streamlines:
        # Overlay streamlines
        lw = speed /0.01 #/speed.max()*0.5
        ax.streamplot(X, Y, mux, muy, density=1.0, linewidth=lw,
                      color='#3399ff', arrowsize=0.1, arrowstyle='->', broken_streamlines=False, zorder=0)

    # Layout bookkeeping
    if plot_title:
        ax.set_title(('vorticity'+' at '+r'$t=%i \Delta t$') % (fr*dt))

    # Set axis limits
    ax.set_aspect('equal')
    ax.set_xlim(0, nx)
    ax.set_ylim(0, ny)
    # Set axis ticks
    if plot_ticks:
        ax.set_xticks([0, nx//2, nx])
        ax.set_xticklabels([0, r'$L/2$', '$L$'], fontsize=8)
        ax.set_yticks([0, ny//2, ny])
        ax.set_yticklabels([0, r'$\dfrac{H}{2}$', '$H$'], fontsize=8)
    else:
        ax.set_xticks(np.arange(0, nx, 50))
        ax.set_xticklabels([])
        ax.set_yticks(np.arange(0, ny, 50))
        ax.set_yticklabels([])
    ax.minorticks_off()
    ax.tick_params(axis='x', which='both', bottom=True,
                   top=False, labelbottom=True, size=2)
    ax.tick_params(axis='y', which='both', left=True,
                   right=False, labelleft=True, size=2)
    ax.margins(0.5, 0.5)

    if (save_plots or save_results) and outdir != None:
        savedir = create_dirs(outdir, save_plots, save_results)
        if save_plots:
            figname = "/vort_frame_%i.png" % (fr)
            fig.savefig(savedir+figname, transparent=True, bbox_inches='tight',
                        pad_inches=0)
        if save_results:
            figname = "/" + \
                os.path.split(outdir)[-1][:-4] + "_vort_frame_%i.png" % (fr)
            fig.savefig(savedir+figname, transparent=True, bbox_inches='tight',
                        pad_inches=0.005)
        plt.close()

    return fig, ax, img


def plot_stress(outdir, fr, dt, figsize, dpi, cmap, cmin, cmax, cvals,
                plot_rmap=False, plot_streamlines=False,
                plot_cbar=False, plot_title=False, plot_ticks=True,
                save_plots=False, save_results=False):
    # Load in the physical fields at frame fr
    nx, ny, rho, ux, uy, speed, vort, Fx, Fy, ss = load_file(outdir, fr)
    # Create mesh
    x = np.arange(0.5, nx+0.5, 1)
    y = np.arange(0.5, ny+0.5, 1)
    X, Y = np.meshgrid(x, y)

    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)
    # Fix axis location
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0,
                    chartBox.width,
                    chartBox.height])
    # Plot density
    img = ax.imshow(ss, cmap=cmap, interpolation='none', vmin=cmin, vmax=cmax,
                    origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])
    # Colorbar
    if plot_cbar:
        cbar = fig.colorbar(img, ax=ax, format="%.4f", pad=0.02, ticks=cvals)
        cbar.ax.set_yticklabels(cvals)

    # Load in the solid at frame fr
    nobjs, remapX_dict, remapY_dict, phi_dict, cc_dict, detF_dict = load_mmap_file(
        outdir, fr)
    # Reference map contour interval
    xlevels = np.arange(0, nx, 5)
    ylevels = np.arange(0, ny, 5)
    # Plot individual solid
    for n in range(nobjs):
        # Plot solid–fluid interface
        ax.contour(X, Y, phi_dict[n], levels=[
                   0.], colors='white', linewidths=0.5, linestyles='solid')
        if plot_rmap:
            # Create reference map mask
            remap_mask = (phi_dict[n] > 0.).astype(int)
            mask = (remap_mask).astype(bool)
            remapX = np.ma.array(remapX_dict[n], mask=mask)
            remapY = np.ma.array(remapY_dict[n], mask=mask)
            # Plot reference map contours
            ax.contour(X, Y, remapX, levels=xlevels, colors='k',
                       linewidths=0.25, linestyles='solid', corner_mask=True)
            ax.contour(X, Y, remapY, levels=ylevels, colors='k',
                       linewidths=0.25, linestyles='solid', corner_mask=True)

    # Layout bookkeeping
    if plot_title:
        ax.set_title(('stress'+' at '+r'$t=%i \Delta t$') % (fr*dt))

    # Set axis limits
    ax.set_aspect('equal')
    ax.set_xlim(0, nx)
    ax.set_ylim(0, ny)
    # Set axis ticks
    if plot_ticks:
        ax.set_xticks([0, nx//2, nx])
        ax.set_xticklabels([0, r'$L/2$', '$L$'], fontsize=8)
        ax.set_yticks([0, ny//2, ny])
        ax.set_yticklabels([0, r'$\dfrac{H}{2}$', '$H$'], fontsize=8)
    else:
        ax.set_xticks([0, nx//2, nx])
        ax.set_xticklabels([])
        ax.set_yticks([0, ny//2, ny])
        ax.set_yticklabels([])
    ax.tick_params(axis='x', which='both', bottom=True,
                   top=False, labelbottom=True)
    ax.tick_params(axis='y', which='both', left=True,
                   right=False, labelleft=True)
    ax.margins(0.5, 0.5)

    if (save_plots or save_results) and outdir != None:
        savedir = create_dirs(outdir, save_plots, save_results)
        if save_plots:
            figname = "/stress_frame_%i.png" % (fr)
            fig.savefig(savedir+figname, transparent=True, bbox_inches='tight',
                        pad_inches=0)
        if save_results:
            figname = "/" + \
                os.path.split(outdir)[-1][:-4] + "_stress_frame_%i.png" % (fr)
            fig.savefig(savedir+figname, transparent=True, bbox_inches='tight',
                        pad_inches=0.005)
        plt.close()

    return fig, ax, img


def plot_trajectory(outdir, fr_end, skip, figsize, dpi, cmap, cmin, cmax, cvals,
                    cm_red, cm_yellow,
                    plot_rmap=False, plot_streamlines=False,
                    plot_cbar=False, plot_title=False, plot_ticks=True,
                    save_plots=False, save_results=False):
    # Create mesh
    nx, ny, rho, ux, uy, speed, vort, Fx, Fy, ss = load_file(outdir, fr_end)
    x = np.arange(0.5, nx+0.5, 1)
    y = np.arange(0.5, ny+0.5, 1)
    X, Y = np.meshgrid(x, y)

    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)
    # Fix axis location
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0,
                    chartBox.width,
                    chartBox.height])

    # Normalize solid stress values
    nss = ss / np.max(ss)

    # Plot solid stress
    img = ax.imshow(nss, cmap=cmap, interpolation='none', vmin=cmin, vmax=cmax,
                    origin='lower', extent=[x.min(), x.max(), y.min(), y.max()], zorder=0)
    # Colorbar
    if plot_cbar:
        cbar = fig.colorbar(img, ax=ax, format="%.4f", pad=0.02, ticks=cvals)
        cbar.ax.set_yticklabels(cvals)

    # Reference map contour interval
    xlevels = np.arange(0, nx, 5)
    ylevels = np.arange(0, ny, 5)

    # Define color
    colors_red = cm_red(np.linspace(0, 1, fr_end+1))
    colors_yellow = cm_yellow(np.linspace(0, 1, fr_end+1))

    # Plot individual solid evolution
    for fr in range(0, fr_end+1, skip):
        # Load in the solid at frame fr
        nobjs, remapX_dict, remapY_dict, phi_dict = load_mmap_file(
            outdir, fr)[:4]

        # TODO: Make the solid color generalizable
        rgba_color = colors_red[fr, :]
        if np.min(rho) < 1.0 and np.max(rho) < 1.125:
            rgba_color = colors_yellow[fr]
        color = matplotlib.colors.to_hex(rgba_color, keep_alpha=True)
        lw = (0.5/fr_end)*fr
        # Loop through solids
        for n in range(nobjs):
            # Plot solid–fluid interface
            ax.contour(X, Y, phi_dict[n], levels=[0.],
                       colors=color, linewidths=lw, linestyles='solid', zorder=1)

    # Plot reference map at final frame
    # Create reference map mask
    remap_mask = (phi_dict[n] > 0.).astype(int)
    mask = (remap_mask).astype(bool)
    remapX = np.ma.array(remapX_dict[n], mask=mask)
    remapY = np.ma.array(remapY_dict[n], mask=mask)
    # Plot masked solid stress
    mss = np.ma.array(nss, mask=mask)
    ax.imshow(mss, cmap=cmap, interpolation='none', vmin=cmin, vmax=cmax,
              origin='lower', extent=[x.min(), x.max(), y.min(), y.max()], zorder=100)

    # Plot reference map contours
    ax.contour(X, Y, phi_dict[n], levels=[0.],
               colors=color, linewidths=0.5, linestyles='solid', zorder=100)
    ax.contour(X, Y, remapX, levels=xlevels, colors=color,
               linewidths=0.25, linestyles='solid', corner_mask=True, zorder=100)
    ax.contour(X, Y, remapY, levels=ylevels, colors=color,
               linewidths=0.25, linestyles='solid', corner_mask=True, zorder=100)

    # Set axis limits
    ax.set_aspect('equal')
    ax.set_xlim(0, nx)
    ax.set_ylim(0, ny)
    # Set axis ticks
    if plot_ticks:
        ax.set_xticks([0, nx//2, nx])
        ax.set_xticklabels([0, r'$L/2$', '$L$'], fontsize=8)
        ax.set_yticks([0, ny//2, ny])
        ax.set_yticklabels([0, r'$\dfrac{H}{2}$', '$H$'], fontsize=8)
    else:
        ax.set_xticks([0, nx//2, nx])
        ax.set_xticklabels([])
        ax.set_yticks([0, ny//2, ny])
        ax.set_yticklabels([])
    ax.tick_params(axis='x', which='both', bottom=True,
                   top=False, labelbottom=True)
    ax.tick_params(axis='y', which='both', left=True,
                   right=False, labelleft=True)
    ax.margins(0.5, 0.5)

    if (save_plots or save_results) and outdir != None:
        savedir = create_dirs(outdir, save_plots, save_results)
        if save_plots:
            figname = "/solid_trajectory.png"
            fig.savefig(savedir+figname, transparent=True, bbox_inches='tight',
                        pad_inches=0)
        if save_results:
            figname = "/" + \
                os.path.split(outdir)[-1][:-4] + "_solid_trajectory.png"
            fig.savefig(savedir+figname, transparent=True, bbox_inches='tight',
                        pad_inches=0.005)
        plt.close()

    return fig, ax, img


def compute_avg_locs_vels(outdir, fr, eps):
    # Load velocity field
    nx, ny, rho, ux, uy, speed = load_file(outdir, fr)[:6]
    # Load solids
    nobjs, remapX_dict, remapY_dict, phi_dict, cc_dict, detF_dict = load_mmap_file(
        outdir, fr)

    # Create mesh
    x = np.arange(0.5, nx+0.5, 1)
    y = np.arange(0.5, ny+0.5, 1)
    X, Y = np.meshgrid(x, y)

    # Compute centroids and weighted average velocity of solids
    avg_locs, avg_vels = [], []
    for n in range(nobjs):
        # Compute Heaviside
        maskphi = heaviside(np.ma.array(phi_dict[n]), eps)
        # Create reference map mask
        remap_mask = (phi_dict[n] > eps).astype(int)
        mask = (remap_mask).astype(bool)
        remapX = np.ma.array(X, mask=mask)
        remapY = np.ma.array(Y, mask=mask)
        # Create velocity mask
        speed_mask = (phi_dict[n] < 0.).astype(int)
        mask = (speed_mask).astype(bool)
        mux = np.ma.array(ux, mask=mask)
        muy = np.ma.array(uy, mask=mask)
        # Compute weighted average position
        cpx = np.sum(np.multiply(remapX, maskphi))/np.sum(maskphi)
        cpy = np.sum(np.multiply(remapY, maskphi))/np.sum(maskphi)
        # Compute weighted average velocity
        cux = np.sum(np.multiply(mux, maskphi))/np.sum(maskphi)
        cuy = np.sum(np.multiply(muy, maskphi))/np.sum(maskphi)
        # Store centroid position
        avg_locs.append([cpx, cpy])
        avg_vels.append([cux, cuy])

    return np.array(avg_locs), np.array(avg_vels)


def heaviside(phi, eps):
    for j in range(phi.shape[0]):
        for i in range(phi.shape[1]):
            if phi[j, i] >= eps:
                phi[j, i] = 0
            elif phi[j, i] <= -eps:
                phi[j, i] = 1
            else:
                phi[j, i] = 1.0-0.5*(1+phi[j, i]/eps+1. /
                                     np.pi*np.sin(np.pi*phi[j, i]/eps))
    return phi
