import matplotlib
import pylab
import matplotlib.patheffects as PathEffects
import matplotlib.transforms as transforms

def annotate(regions, lobediv, L, ax, zorder=0):
    transY = transforms.blended_transform_factory(ax.transAxes, ax.transData)
    transX = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    clr = (0.4,0.4,0.4)
    for n,(region, start, end) in enumerate(regions):
        ax.fill([start+0.5, start+0.5, end+0.5, end+0.5], 
                [+0.5,L+0.5,L+0.5,+0.5],
                fill=False, hatch='\\' , ec=clr, lw=0, zorder=zorder+1)
        ax.fill([0.5,L+0.5,L+0.5,0.5],
                [start+0.5, start+0.5, end+0.5, end+0.5], 
                fill=False, hatch='/', ec=clr, lw=0, zorder=zorder+1)
        ax.axvline(start+0.5, color=clr, zorder=zorder+1).set_snap(False)
        ax.axvline(end+0.5, color=clr, zorder=zorder+1).set_snap(False)
        ax.axhline(start+0.5, color=clr, zorder=zorder+1).set_snap(False)
        ax.axhline(end+0.5, color=clr, zorder=zorder+1).set_snap(False)
        effects = PathEffects.withStroke(linewidth=3, foreground="w")
        ax.text(0.990, end+0.5, region, 
                verticalalignment='top', horizontalalignment='right', 
                transform=transY, zorder=100, path_effects=[effects])
        if region != 'Activation Loop':
            ax.text(end, 0.990, region, rotation='vertical', 
                    verticalalignment='top', horizontalalignment='right', 
                    transform=transX, zorder=100, path_effects=[effects])
    
    ax.plot([+0.5, L+0.5], [+0.5, L+0.5], 'k-', zorder=zorder+1)
    ax.axvline(lobediv+0.5, color='k', zorder=zorder+1).set_snap(False)
    ax.axhline(lobediv+0.5, color='k', zorder=zorder+1).set_snap(False)
    ax.set_xlim(+0.5,L+0.5)
    ax.set_ylim(+0.5,L+0.5)
