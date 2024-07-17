include("../LinStab/mjcl_stab.jl")
using .LinStab

# loading modules
using PyPlot, CSV, DataFrames, Printf, FFTW, PlotUtils

using JLD2, FileIO, Statistics, LaTeXStrings

# controls ratio of interfaces' stratification
gammas = round.(collect(range(0.2,3,15)),sigdigits=2) #

# controls ratio of interface shears
alphas = round.(collect(range(1,5,21)),sigdigits=2)

# topo parameters
h0s_ext = collect(range(-10^-3,10^-3,30));
h0s_ext = [h0s_ext[1:15]; 0.; h0s_ext[16:end]];


# h0s_ext = [-10 .^reverse(collect(range(-7,-3,6))); 0. ; 10 .^collect(range(-7,-3,6))]
# h0s_fill = [ - collect(range(10^-4.5,10^-3.2,10)); collect(range(10^-4.5,10^-3.2,10))]
        
# h0s_ext = sort([h0s_ext; h0s_fill]);

h0_x = [-10 .^reverse(collect(range(-7,-3,17))); [0]; 10 .^collect(range(-7,-3,17))]

kts = [0]  # zero wavenumber for sloping topography

Ng = length(gammas); Na = length(alphas); Nh_ext = length(h0s_ext); Nk = length(kts);
        
H = [500.,1000.,2500.]; Nz = length(H);
H_grid = -[0, H[1], sum(H[1:2]), sum(H)];

# get file names
dat_path = "../data/data_jld2_02/"
all_files = readdir(dat_path)

Lx = 1000.e3

filename(gam,alp,h0,k) = "threelayer_power_iter_gamma"*string(gam)*"_alpha"*string(alp)*"_h0"*string(round(h0*Lx,digits=9))*"_kt"*string(Int(k))*"_res256.jld2"

nlayers = 3                 # number of layers
f0 = f₀= 8.3e-5  # 1.0e-4            # Coriolis param [s^-1]
g = 9.81                    # gravity
H = [500., 1000., 2500.]     # the rest depths of each layer; [250., 500., 3000.]

H32 = 0.5 * (H[1] + H[2])
H52 = 0.5 * (H[2] + H[3])

U = [0.05,0.025,0.0]

# setting base density profile
rho = ρ = [0.0, 1025.0, 1025.75]         # the density of each layer

U1 = U[2] .+ alphas*(H32/H52)*(U[2]-U[3])

rho1g3 = rho[2] - gammas[3] * (rho[3] - rho[2]) * (H32/H52)

rho1g13 = rho[2] - gammas[13] * (rho[3] - rho[2]) * (H32/H52)

S52g3 = -(f0*rho1g3/g)*((U[2] .- U[3])./(rho[2]-rho[3]))

S52g13 =  -(f0*rho1g13/g)*((U[2] .- U[3])./(rho[2]-rho[3]))

alt_x = (H[2]/H[3])*((g*(rho1g3 - rho[3])./rho1g3) ./ (f0*(U[2] .- U1))) * h0s_ext' # Wright (1980) eq. at beginning of section 4b

dat_in = load(dat_path*filename(gammas[3],alphas[1],h0s_ext[1],kts[1]));
k_xr = dat_in["csv_data"]["k_growth_lsa"][129:end]

Nk = length(k_xr)
Lx = 1000.e3

############################################################################
function alpha_func(h0,h0s)
    eps = 1.e-9
    s1 = abs.(log10.(abs.(h0s))) .- minimum(abs.(log10.(abs.(h0s))))
    alphav = (abs(log10(abs(h0+eps))) - minimum(abs.(log10.(abs.(h0s)))) ) / maximum(s1)
    return (1 - alphav)^3
end
    
alpha_func(h0s_ext[14],h0s_ext)

using PlotUtils

cm = cgrad(:coolwarm);
colors = [cm[i] for i in range(0,1,100)]

hc_range = collect(range(h0s_ext[1],h0s_ext[end],100))

function h_color(h_in,hc_range,colors)
    c_ind = argmin(abs.(h_in .- (hc_range)))

    color_out = [red(colors[c_ind]),green(colors[c_ind]),blue(colors[c_ind])]
    return color_out
end

fig,ax = PyPlot.subplots(1,1,figsize=(9,7.5))
fig.tight_layout(pad=0.0)
ct = 1

# plot k_growth_emp as function of h0 and alpha for a given gamma; Charney-favoring
# kspec_gam3 = zeros(Nh_ext,Na)

gam = gammas[13];
kt = kts[1];
alp = alphas[19];

for (i,h0) in enumerate(h0s_ext)
#     for (j,alp) in enumerate(alphas)
        fil = findall(x->x==filename(gam,alp,h0,kt),all_files)
        if isempty(fil)==false
            dat_in = load(dat_path*filename(gam,alp,h0,kt));
            a = dat_in["csv_data"]["psi2_full"]
            b = mean(abs.(rfft(a,[2])),dims=1)' ;

            if h0<0
                ax.plot(b,"-",color=h_color(h0,hc_range,colors),alpha=alpha_func(h0,h0s_ext))
            else
                ax.plot(b,color=h_color(h0,hc_range,colors),alpha=alpha_func(h0,h0s_ext))
            end

            color=h_color(h0,hc_range,colors)
        end
#     end
end

ax.set_xlim(0.001,20)


#############################################################################
function plot_spec(rd,gam,alp,kt,ax1,ax2,xmax,ymax,l_flag)
    for (i,h0) in enumerate(h0s_ext)
    #     for (j,alp) in enumerate(alphas)
            fil = findall(x->x==filename(gam,alp,h0,kt),all_files)
            if isempty(fil)==false
                dat_in = load(dat_path*filename(gam,alp,h0,kt));
                b = dat_in["csv_data"]["CV32"][end][1][1:end-1]
                c = dat_in["csv_data"]["CV52"][end][1][1:end-1]
                if h0<0.
                    ax1.plot(k_xr*rd,b,"-",color=h_color(h0,hc_range,colors),linewidth=2) # ,alpha=alpha_func(h0,h0s_ext))
                    ax2.plot(k_xr*rd,c,"-",color=h_color(h0,hc_range,colors),linewidth=2) # ,alpha=alpha_func(h0,h0s_ext))
#                 elseif i == 16 || i == 15
#                     ax1.plot(k_xr*rd,b,color="black",linewidth=2) # ,alpha=alpha_func(h0,h0s_ext))
#                     ax2.plot(k_xr*rd,c,color="black",linewidth=2) # ,alpha=alpha_func(h0,h0s_ext))
                elseif h0==0.
#                     ax1.plot(k_xr*rd,b,color="black",linewidth=2) # ,alpha=alpha_func(h0,h0s_ext))
#                     ax2.plot(k_xr*rd,c,color="black",linewidth=2) # ,alpha=alpha_func(h0,h0s_ext))
                else
                    ax1.plot(k_xr*rd,b,color=h_color(h0,hc_range,colors),linewidth=2) # ,alpha=alpha_func(h0,h0s_ext))
                    ax2.plot(k_xr*rd,c,color=h_color(h0,hc_range,colors),linewidth=2) # ,alpha=alpha_func(h0,h0s_ext))
                end

                color=h_color(h0,hc_range,colors)
            end
    #     end
    end

    zero_ind = argmin(abs.(h0s_ext))
    dat_in = load(dat_path*filename(gam,alp,h0s_ext[zero_ind],kt));
    b = dat_in["csv_data"]["CV32"][end][1][1:end-1]
    c = dat_in["csv_data"]["CV52"][end][1][1:end-1]
    ax1.plot(k_xr*rd,b,color="black",linewidth=2) # ,alpha=alpha_func(h0,h0s_ext))
    ax2.plot(k_xr*rd,c,color="black",linewidth=2) # ,alpha=alpha_func(h0,h0s_ext))

    ax1.text(0.4, 0.65, L"\frac{S_{3/2}}{S_{5/2}} = "*latexstring(round(alp/gam,sigdigits=2)), fontsize=26, horizontalalignment="center",verticalalignment="center", transform = ax1.transAxes)

    ax1.set_xlim(0.001,xmax)
    ax2.set_xlim(0.001,xmax)

    ax1.set_ylim(0,ymax)
    ax2.set_ylim(0,ymax)

    ax1.set_xticklabels([])

    if l_flag==0
        ax1.set_yticklabels([])
        ax2.set_yticklabels([])
    else
        ax1.ticklabel_format(axis="y", scilimits=(0,0))
        ax2.ticklabel_format(axis="y", scilimits=(0,0))

        ax1.yaxis.offsetText.set(size=18)
        ax2.yaxis.offsetText.set(size=18)
    end

    ax1.tick_params(labelsize=18)
    ax2.tick_params(labelsize=18)

#     ax1.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=true))
#     ax2.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=true))
    ax1.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.5))
    ax2.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.5))

end

#################################################
# manually adding deformation radii for plotting
rd26 = 34049
rd06 = 31744



#################################################
matplotlib[:rc]("text", usetex= "true")
matplotlib[:rc]("text.latex", preamble= "\\usepackage{amsmath,amsthm}")


using PlotUtils

f_size = 26

cm = cgrad(:coolwarm);
colors = [cm[i] for i in range(0,1,100)]

hc_range = collect(range(h0s_ext[1],h0s_ext[end],100))

function h_color(h_in,hc_range,colors)
    c_ind = argmin(abs.(h_in .- (hc_range)))

    color_out = [red(colors[c_ind]),green(colors[c_ind]),blue(colors[c_ind])]
    return color_out
end

# fig,ax = PyPlot.subplots(2,3,figsize=(10,9))
fig,ax = PyPlot.subplots(2,4,figsize=(12,7))
fig.tight_layout(pad=0.0)
# ax1=ax[0][0]; ax2=ax[1][0]; ax3=ax[0][1]; ax4=ax[1][1]; ax5=ax[0][2]; ax6=ax[1][2]; ax7=ax[0][3]; ax8=ax[1][3];
ax1=ax[1]; ax2=ax[2]; ax3=ax[3]; ax4=ax[4]; ax5=ax[5]; ax6=ax[6]; ax7=ax[7]; ax8=ax[8];

plot_spec(rd26,gammas[13],alphas[3],kts[1],ax1,ax2,2.25,0.025,1)
plot_spec(rd26,gammas[13],alphas[9],kts[1],ax3,ax4,2.25,0.025,0)
plot_spec(rd26,gammas[13],alphas[11],kts[1],ax5,ax6,2.25,0.025,0)
plot_spec(rd06,gammas[3],alphas[1],kts[1],ax7,ax8,2.25,0.025,0)

# plot_spec(gammas[3],alphas[1],kts[1],ax3,ax4,18,0.025)
# plot_spec(gammas[13],alphas[end],kts[1],ax5,ax6,18,0.25)


ax1.set_ylabel(L"\widehat{C}_{\mathrm{V},3/2} \ \ [\mathrm{m}^4 \, \mathrm{s}^{-3}]",fontsize=f_size)
ax2.set_ylabel(L"\widehat{C}_{\mathrm{V},5/2} \ \ [\mathrm{m}^4 \, \mathrm{s}^{-3}]",fontsize=f_size)

# ax4.set_xlabel(L"k_x \cdot L_x / ( 2 \pi ) ",fontsize=20)

fig.text(0.375, -0.05, L"k_x L_{\mathrm{d,} 1} ", va="center", rotation="horizontal", fontsize=f_size)


# ax6.text(-0.25, 0.9, L"\text{linear y axes}", fontsize=26, horizontalalignment="center",verticalalignment="center",color="red", transform = ax6.transAxes);
# ax6.text(-0.25, 0.7, L"\text{all with equal range}", fontsize=26, horizontalalignment="center",verticalalignment="center",color="red", transform = ax6.transAxes);

fig.subplots_adjust(right=0.8)

cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])

cm = matplotlib.cm.coolwarm
norm = matplotlib.colors.Normalize(vmin=h0s_ext[1]/S52g13, vmax=h0s_ext[end]/S52g13)

cb1 = matplotlib.colorbar.ColorbarBase(cbar_ax, cmap=cm,
                                norm=norm,
                                orientation="vertical")
cb1.ax.set_ylabel(L"\frac{\partial_y \eta_b }{ S_{5/2} }",size=36,rotation=360,y=0.6,horizontalalignment="left") # ,y=0.475)

cb1.formatter.set_powerlimits((-3,3))
cb1.ax.tick_params(labelsize=20)
cb1.ax.yaxis.offsetText.set(size=20)
cb1.ax.axhline(y=0., c="black",linewidth=3.)


ax1.text(0.15, 0.9, L"\mathbf{(a)}", fontsize=30, horizontalalignment="center",verticalalignment="center", transform = ax1.transAxes);
ax2.text(0.15, 0.9, L"\mathbf{(e)}", fontsize=30, horizontalalignment="center",verticalalignment="center", transform = ax2.transAxes);
ax3.text(0.15, 0.9, L"\mathbf{(b)}", fontsize=30, horizontalalignment="center",verticalalignment="center", transform = ax3.transAxes);
ax4.text(0.15, 0.9, L"\mathbf{(f)}", fontsize=30, horizontalalignment="center",verticalalignment="center", transform = ax4.transAxes);
ax5.text(0.15, 0.9, L"\mathbf{(c)}", fontsize=30, horizontalalignment="center",verticalalignment="center", transform = ax5.transAxes);
ax6.text(0.15, 0.9, L"\mathbf{(g)}", fontsize=30, horizontalalignment="center",verticalalignment="center", transform = ax6.transAxes);
ax7.text(0.15, 0.9, L"\mathbf{(d)}", fontsize=30, horizontalalignment="center",verticalalignment="center", transform = ax7.transAxes);
ax8.text(0.15, 0.9, L"\mathbf{(h)}", fontsize=30, horizontalalignment="center",verticalalignment="center", transform = ax8.transAxes);

ax1.set_title(L"\mathrm{Bottom-intensified}", fontsize=20,y=1.065)
ax3.set_title(L"\mathrm{Full-column}", fontsize=20,y=1.065)
ax5.set_title(L"\mathrm{Mixed-type}", fontsize=20,y=1.065)
ax7.set_title(L"\mathrm{Surface-intensified}", fontsize=20,y=1.065)


ax1.text(0.7, 0.9, L"\gamma=2.6", color="red", fontsize=22, horizontalalignment="center",verticalalignment="center", transform = ax1.transAxes);
ax3.text(0.7, 0.9, L"\gamma=2.6", color="red", fontsize=22, horizontalalignment="center",verticalalignment="center", transform = ax3.transAxes);
ax5.text(0.7, 0.9, L"\gamma=2.6", color="red", fontsize=22, horizontalalignment="center",verticalalignment="center", transform = ax5.transAxes);
ax7.text(0.7, 0.9, L"\gamma=0.6", color="red", fontsize=22, horizontalalignment="center",verticalalignment="center", transform = ax7.transAxes);

savefig("./vert_trans_ex.pdf",bbox_inches="tight")















