# plot for flat bottom alpha-gamma space

# load modules
include("../LinStab/mjcl_stab.jl")
using .LinStab

using PyPlot, Measures, LaTeXStrings, JLD2, FileIO, Statistics

matplotlib[:rcParams]["axes.unicode_minus"]=false
PyPlot.matplotlib[:rc]("hatch",linewidth=1.5)
PyPlot.matplotlib[:rc]("hatch",color="red")

PyPlot.matplotlib[:rc]("text", usetex=true) 
PyPlot.matplotlib[:rc]("text.latex", preamble= "\\usepackage{amsmath,amsthm}")
PyPlot.matplotlib[:rc]("figure", autolayout=false)


pygui(true)

# set static vars
L = Lx = Ly = 1000.e3                   # domain size [m]
beta = β = 0 #1.14052e-11              # the y-gradient of planetary PV
Ny = Nx = n = 256               # 2D resolution = n²
eta = 0

nlayers = 3                 # number of layers
f0 = f₀= 8.3e-5            # Coriolis param [s^-1] 
g = 9.81                    # gravity
H = [500., 1000., 2500.]     # the rest depths of each layer; [250., 500., 3000.] 

# H = [1000., 1000., 1000.]
# H = [1000., 500., 2500.]
Nz = length(H)

U0 = [0.05,0.025,0.0]

alphas = collect(range(0.6,5.0,101))

Na = length(alphas)

# set gamma range
gammas = collect(range(0.6,5.0,101))

Ng = length(gammas)

# setting density profile as function of gamma
rho = ρ = [1024.0, 1025.0, 1025.75]         # the density of each layer

V = zeros(nlayers)

H32 = 0.5 * (H[1] + H[2])
H52 = 0.5 * (H[2] + H[3])

##################################################

a = load("./LSA_data/ag_kmax.jld")
max_k = a["max_k"]

a = load("./LSA_data/ag_phi12.jld")
phase_shift_12 = a["phase_shift_12"]

a = load("./LSA_data/ag_phi32.jld")
phase_shift_32 = a["phase_shift_32"]

a = load("./LSA_data/ag_sigma_max.jld")
sigma_gamma_alpha_max = a["sigma_gamma_alpha_max"]

Qy_gamma_alpha  = zeros((Ng,Na,Nz))

rigid_lid=true

for i=range(1,Ng)
    for j=range(1,Na)
        gamma = gammas[i]
        alpha = alphas[j]

        U0[1] = U0[2] + alpha*(H32/H52)*(U0[2]-U0[3])

        rho[1] = ρ[1] = rho[2] - gamma * (rho[3] - rho[2]) * (H32/H52)

        N_xr = round(Int,Nx/2)

        S = LinStab.calc_stretching_mat(Nz,rho,f0,H,rho[1],g,rigid_lid,eta)

        # change dimensions of U and V to match domain size
        U2 = zeros(1,Nz); U2[:] = U0; U2 = repeat(U2,outer=(Ny,1))
        U = zeros(1,Ny,Nz); U[1,:,:] = U2

        a = LinStab.calc_PV_grad_y(U,beta,(f0/H[end])*eta,Ny,Nz,0,S)
        
        Qy_gamma_alpha[i,j,:] = mean(a[1,:,:],dims=1)
    end
end


# building filter
using ImageFiltering
blur_12 = imfilter(phase_shift_12, Kernel.gaussian(1))

##############################################

msk = (Qy_gamma_alpha[:,:,2] .< Qy_gamma_alpha[:,:,3]) 
# msk =  (Qy_gamma_alpha[:,:,1] .> - Qy_gamma_alpha[:,:,2])  # (Qy_gamma_alpha[:,:,2] .> Qy_gamma_alpha[:,:,3]) .&& 
msk = convert(Matrix{Union{Float64, Bool}}, msk)
msk[msk .== 0.0] .= NaN

msk2 = (H[2] .* Qy_gamma_alpha[:,:,2] .< H[3] .* Qy_gamma_alpha[:,:,3]) 
msk2 = convert(Matrix{Union{Float64, Bool}}, msk2)
msk2[msk2 .== 0.0] .= NaN


#######################################################################
# A COMBINED FIGURE 
fig,ax = PyPlot.subplots(2,2,figsize=(14,14))
fig.tight_layout(pad=3.0)
ax2 = ax[1]; ax3 = ax[3];
ax1 = ax[2]; ax4 = ax[4];

f_size = 40.
l_size = 28.

ps_min = minimum([phase_shift_12 ./ (pi) ; phase_shift_32 ./ (pi)])
ps_max = maximum([phase_shift_12 ./ (pi) ; phase_shift_32 ./ (pi)])

pc2=ax2.contour(gammas,alphas,(blur_12 ./ (pi))',linewidths=3.,levels=[0.05,0.08,0.11,0.14,0.17,0.2,0.23,0.26,0.29], cmap=PyPlot.cm.Greys,vmin=-0.1)
ax2.tick_params(labelsize=l_size)
ax2.set_xticklabels([])
ax2.set_ylim([alphas[1]-0.02,alphas[end]+0.02])
ax2.set_xlim([gammas[1]-0.02,gammas[end]+0.02])
ax2.set_title(L"\Delta \phi_{3/2}", fontsize = f_size)
ax2.pcolor(gammas,alphas, msk',cmap="Greens",vmin=0., vmax=8.0)
ax2.set_xticks(ax2.get_yticks()[2:end-1])
ax2.text(0.1, 0.9, L"\mathbf{(a)}", fontsize=50, horizontalalignment="center",verticalalignment="center", transform = ax2.transAxes);
cl2 = ax2.clabel(pc2, inline=true, fontsize=20, manual=[(3,1),(3,1.85),(2.5,2),(2.25,2.9),(2.2,3.0),(4,3.5),(4,4),(2,4),(1.5,4),(1,4)])
# [txt.set_backgroundcolor("white") for txt in cl2]
# [txt.set_bbox(dict(facecolor="white", edgecolor="none", pad=0)) for txt in cl2]


pc3=ax3.contour(gammas,alphas,(phase_shift_32 ./ (pi))',linewidths=3.,levels=[0.02,0.1,0.18,0.28], cmap=PyPlot.cm.Greys, vmin=-0.4, vmax=0.34) # [0.05,0.1,0.15,0.2,0.25,0.3,0.35]
ax3.tick_params(labelsize=l_size)
ax3.set_ylim([alphas[1]-0.02,alphas[end]+0.02])
ax3.set_xlim([gammas[1]-0.02,gammas[end]+0.02])
ax3.set_yticklabels([])
ax3.set_xticklabels([])
ax3.set_title(L"\Delta \phi_{5/2}", fontsize = f_size)
ax3.set_xticks(ax2.get_yticks()[2:end-1])
ax3.pcolor(gammas,alphas, msk', cmap="Greens",vmin=0., vmax=8.0)
ax3.text(0.1, 0.9, L"\mathbf{(b)}", fontsize=50, horizontalalignment="center",verticalalignment="center", transform = ax3.transAxes);
cl3 = ax3.clabel(pc3, pc3.levels, inline=true, fontsize=20, manual=[(1.75,4),(2.75,4.0),(3.5,4.25),(4.25,3),(4.25,2)])
# [txt.set_backgroundcolor("white") for txt in cl3]
# [txt.set_bbox(dict(facecolor="white", edgecolor="none", pad=0)) for txt in cl3]

pc1=ax1.contour(gammas,alphas,(2*pi)^-1 *Lx*max_k',colors="black",linewidths=3.,levels=[2,3,4,5,6,7,8])
ax1.pcolor(gammas,alphas, msk',cmap="Greens",vmin=0., vmax=8.0)
ax1.tick_params(labelsize=l_size)
ax1.set_title(L"k_{x\mathrm{, max}} L_x / (2 \pi)", fontsize =f_size)
ax1.set_ylim([alphas[1],alphas[end]])
ax1.set_xlim([gammas[1],gammas[end]])
ax1.text(0.1, 0.9, L"\mathbf{(c)}", fontsize=50, horizontalalignment="center",verticalalignment="center", transform = ax1.transAxes);
ax1.set_xticks(ax1.get_yticks()[2:end])
cl1 = ax1.clabel(pc1, pc1.levels, inline=true, fontsize=20, manual=[(4,2.5),(2.9,4.5),(2.15,3.75),(1.25,3.75),(1.1,3),(0.75,2.25)]) # for 4 (2.25,4.2)
# [txt.set_backgroundcolor("white") for txt in cl1]
# [txt.set_bbox(dict(facecolor="white", edgecolor="none", pad=0)) for txt in cl1]


pc4=ax4.contour(gammas,alphas,log10.(sigma_gamma_alpha_max)',colors="black",linewidths=3.,levels=[-6.8,-6.725,-6.65,-6.575,-6.5,-6.425,-6.35,-6.275,-6.2,-6.125,-6.05])
ax4.tick_params(labelsize=l_size)
ax4.pcolor(gammas,alphas, msk',cmap="Greens",vmin=0., vmax=8.0)
ax4.set_yticklabels([])
ax4.set_ylim([alphas[1],alphas[end]])
ax4.set_xlim([gammas[1],gammas[end]])
ax4.set_title(L"\mathrm{log}_{10} (\sigma_\mathrm{max} \cdot 1 \mathrm{s})", fontsize = f_size)
ax4.text(0.1, 0.9, L"\mathbf{(d)}", fontsize=50, horizontalalignment="center",verticalalignment="center", transform = ax4.transAxes);
ax4.set_xticks(ax1.get_yticks()[2:end])
cl4 = ax4.clabel(pc4, pc4.levels[1:2:end], inline=true, fontsize=20)
# [txt.set_backgroundcolor("white") for txt in cl4]
# [txt.set_bbox(dict(facecolor="white", edgecolor="none", pad=0)) for txt in cl4]


fig.text(0.5, 0., L"\gamma", ha="center",fontsize=f_size)
fig.text(0., 0.5, L"\alpha", ha="center",fontsize=f_size, rotation="vertical")

# savefig("./figs/flat_comb_cont.png",bbox_inches="tight")



#######################################################################
# A COMBINED FIGURE 
fig,ax = PyPlot.subplots(2,2,figsize=(7.5,6))
fig.tight_layout(pad=1.0)
ax2 = ax[1]; ax3 = ax[3];
ax1 = ax[2]; ax4 = ax[4];

f_size = 22.
l_size = 14.

ps_min = minimum([phase_shift_12 ./ (pi) ; phase_shift_32 ./ (pi)])
ps_max = maximum([phase_shift_12 ./ (pi) ; phase_shift_32 ./ (pi)])

pc2=ax2.pcolormesh(gammas,alphas,(phase_shift_12 ./ (pi))',vmin=ps_min,vmax=ps_max,cmap=PyPlot.cm.Blues,linewidth=0,rasterized=true)
ax2.tick_params(labelsize=l_size)
ax2.set_xticklabels([])
ax2.set_ylim([alphas[1]-0.02,alphas[end]+0.02])
ax2.set_xlim([gammas[1]-0.02,gammas[end]+0.02])
ax2.set_title(L"\Delta \phi_{3/2}", fontsize = f_size, y=1.025)
ax2.pcolor(gammas,alphas, msk', hatch="/", alpha=0.)
ax2.set_xticks(ax2.get_yticks()[2:end-1])
ax2.text(0.1, 0.9, L"\mathbf{(a)}", fontsize=25, horizontalalignment="center",verticalalignment="center", transform = ax2.transAxes);

cb = fig.colorbar(pc2)
cb.ax.tick_params(labelsize=l_size)
cb.ax.yaxis.offsetText.set(size=l_size)



pc3=ax3.pcolormesh(gammas,alphas,(phase_shift_32 ./ (pi))',vmin=ps_min,vmax=ps_max,cmap=PyPlot.cm.Blues,linewidth=0,rasterized=true)
ax3.tick_params(labelsize=l_size)
ax3.set_ylim([alphas[1]-0.02,alphas[end]+0.02])
ax3.set_xlim([gammas[1]-0.02,gammas[end]+0.02])
ax3.set_yticklabels([])
ax3.set_xticklabels([])
ax3.set_title(L"\Delta \phi_{5/2}", fontsize = f_size, y=1.025)

ax3.set_xticks(ax2.get_yticks()[2:end-1])
ax3.pcolor(gammas,alphas, msk', hatch="/", alpha=0.)
ax3.text(0.1, 0.9, L"\mathbf{(b)}", fontsize=25, horizontalalignment="center",verticalalignment="center", transform = ax3.transAxes);


cb = fig.colorbar(pc3)
cb.ax.tick_params(labelsize=l_size)
cb.ax.yaxis.offsetText.set(size=l_size)

pc1=ax1.pcolormesh(gammas,alphas,(2*pi)^-1 *Lx*max_k',cmap=PyPlot.cm.summer,alpha=1.0,linewidth=0,rasterized=true)   # cmap=PyPlot.cm.cividis
ax1.pcolor(gammas,alphas, msk', hatch="/", alpha=0.)
ax1.tick_params(labelsize=l_size)
ax1.set_title(L"k_{x\mathrm{, max}} L_x / (2 \pi)", fontsize =f_size, y=1.025)
ax1.set_ylim([alphas[1],alphas[end]])
ax1.set_xlim([gammas[1],gammas[end]])

cb = fig.colorbar(pc1)
cb.ax.tick_params(labelsize=l_size)
cb.ax.yaxis.offsetText.set(size=l_size)
ax1.text(0.1, 0.9, L"\mathbf{(c)}", fontsize=25, horizontalalignment="center",verticalalignment="center", transform = ax1.transAxes);
ax1.set_xticks(ax1.get_yticks()[2:end])


pc4=ax4.pcolormesh(gammas,alphas,(sigma_gamma_alpha_max .* 3600 .* 24 .* 7)',cmap=PyPlot.cm.cividis,linewidth=0,rasterized=true) # cmap=PyPlot.cm.plasma; log10.(sigma_gamma_alpha_max)'
ax4.tick_params(labelsize=l_size)
ax4.pcolor(gammas,alphas, msk', hatch="/", alpha=0.)
ax4.set_yticklabels([])
ax4.set_ylim([alphas[1],alphas[end]])
ax4.set_xlim([gammas[1],gammas[end]])
ax4.set_title(L"\sigma_\mathrm{max} \ \ { \scriptstyle  \left [ \mathrm{week}^{-1} \right ] }", fontsize = f_size, y=1.025)
cb = fig.colorbar(pc4)
cb.ax.tick_params(labelsize=l_size)
cb.ax.yaxis.offsetText.set(size=l_size)
ax4.text(0.1, 0.9, L"\mathbf{(d)}", fontsize=25, horizontalalignment="center",verticalalignment="center", transform = ax4.transAxes); # color="#F5F5DC" # pretty off-white
ax4.set_xticks(ax1.get_yticks()[2:end])

fig.text(0.5, 0., L"\gamma", ha="center",fontsize=f_size)
fig.text(0., 0.5, L"\alpha", ha="center",fontsize=f_size, rotation="vertical")

savefig("./figs/flat_comb.pdf",bbox_inches="tight")


savefig("./figs/flat_comb.png",bbox_inches="tight")

