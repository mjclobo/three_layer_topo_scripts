# plot in alpha-eta_b space

# load modules
include("../LinStab/mjcl_stab.jl")
using .LinStab

using PyPlot, Measures, LaTeXStrings, FileIO, JLD2

matplotlib[:rcParams]["axes.unicode_minus"]=false
PyPlot.matplotlib[:rc]("text", usetex=true) 
PyPlot.matplotlib[:rc]("text.latex", preamble= "\\usepackage{amsmath,amsthm}")

pygui(true)

# set static vars
global L = global Lx = global Ly = 1000.e3                   # domain size [m]
beta = β = 0 # 1.14052e-11              # the y-gradient of planetary PV
Ny = Nx = n = 256                # 2D resolution = n²

nlayers = 3                 # number of layers
f0 = f₀= 8.3e-5            # Coriolis param [s^-1] 
g = 9.81                    # gravity
H = [500., 1000., 2500.]     # the rest depths of each layer; [250., 500., 3000.] 
Nz = length(H)

U = [0.05,0.025,0.0]

# controls ratio of interface shears
alphas = round.(collect(range(1,5,101)),sigdigits=2) #

Na = length(alphas)

# set gamma range
gammas = round.(collect(range(0.2,3,15)),sigdigits=2) #
# gammas = round.(collect(range(0.2,3,51)),sigdigits=2) 
Ng = length(gammas)

gamma = 2.6 # gammas[13]  # gammas[4] shows beta destab best; gammas[8] shows instability superposition

# setting density profile as function of gamma
rho = ρ = [1024.0, 1025.0, 1025.75]         # the density of each layer

V = zeros(nlayers)

# Linear bottom drag
μ = 0 # f0*10/2/H[3]            # d_E is a typical value of 5-10 m in the ocean [Weatherly and Martin, 1978; Perlin et al., 2007]

# step thru gamma while storing vertical structure of most unstable mode
Nk = round(Int,Nx/2)
sigma_gamma = zeros((Ng,Nk))

h0s = collect(range(-10^-3,10^-3,101))
Nh = length(h0s)

H32 = 0.5 * (H[1] + H[2])
H52 = 0.5 * (H[2] + H[3])


rho[1] = ρ[1] = rho[2] - gamma * (rho[3] - rho[2]) * (H32/H52)

U1 = U[2] .+ alphas*(H32/H52)*(U[2]-U[3])
rho1g3 = rho[2] - gammas[3] * (rho[3] - rho[2]) * (H32/H52)

rho1g13 = rho[2] - gammas[13] * (rho[3] - rho[2]) * (H32/H52)

S52g3 = -(f0*rho1g3/g)*((U[2] .- U[3])./(rho[2]-rho[3]))

S52g13 =  -(f0*rho1g13/g)*((U[2] .- U[3])./(rho[2]-rho[3]))


#############################################################
# loading data
a = load("./LSA_data/ae_phi12.jld")
phase_shift_12 = a["phase_shift_12"]

a = load("./LSA_data/ae_phi32.jld")
phase_shift_32 = a["phase_shift_32"]

a = load("./LSA_data/ae_Qy_alpha_eta.jld")
Qy_alpha_eta = a["Qy_alpha_eta"]


## filtering
using ImageFiltering
blur_12 = imfilter(phase_shift_12, Kernel.gaussian(1))
blur_32 = imfilter(phase_shift_32, Kernel.gaussian(1))


#######################################################################
# making mask
msk = (Qy_alpha_eta[:,:,2] .< Qy_alpha_eta[:,:,3]) .&& (Qy_alpha_eta[:,:,2] .< 0)
msk = convert(Matrix{Union{Float64, Bool}}, msk)
msk[msk .== 0] .= NaN


#########################################################################
# contours per Isaac's recommendation
fig,ax = PyPlot.subplots(1,2,figsize=(16,7))
fig.tight_layout(pad=5.0)
ax2 = ax[1]; ax3 = ax[2];

f_size = 56.
l_size = 28.

ps_min = 0. #  minimum([phase_shift_12 ./ (2*pi) ; phase_shift_32 ./ (2*pi)])
ps_max = maximum([phase_shift_12 ./ (pi) ; phase_shift_32 ./ (pi)])


pc2=ax2.contour(h0s/S52g13,alphas/gamma,(blur_12 ./ (pi))', levels=[0.02,0.05,0.1,0.15,0.2,0.25], linewidths=3.0, cmap=PyPlot.cm.Greys,vmin=-0.2)  # [0.0,0.06,0.11,0.16,0.21,0.26]
ax2.clabel(pc2, pc2.levels, inline=true, fontsize=20.)
ax2.tick_params(labelsize=l_size)
ax2.set_ylim([alphas[1]/gamma, alphas[end]/gamma])
ax2.set_xlim([h0s[1]/S52g13, h0s[end]/S52g13])
# ax2.set_yticklabels([])
ax2.set_title(L"\Delta \phi_{3/2}", fontsize = f_size)
# ax2.set_xticklabels([])
ax2.plot([h0s[1]/S52g13, h0s[end]/S52g13],[1,1],color="cyan",linestyle="dashed",linewidth=4.)
ax2.pcolor(h0s/S52g13,alphas/gamma, msk', cmap="Greens",vmin=0., vmax=8.0)
ax2.text(0.1, 0.9, L"\mathbf{(a)}", fontsize=50, horizontalalignment="center",verticalalignment="center", transform = ax2.transAxes);



pc3=ax3.contour(h0s/S52g13,alphas/gamma,(blur_32 ./ (pi))', levels=[0.0,0.08,0.16,0.24,0.32,0.4], linewidths=3.0, cmap=PyPlot.cm.Greys,vmin=-0.2) # [0.0,0.075,0.15,0.225,0.3,0.375]
ax3.clabel(pc3, pc3.levels, inline=true, fontsize=20.)
ax3.tick_params(labelsize=l_size)
ax3.set_ylim([alphas[1]/gamma, alphas[end]/gamma])
ax3.set_xlim([h0s[1]/S52g13, h0s[end]/S52g13])
ax3.set_yticklabels([])
ax3.set_title(L"\Delta \phi_{5/2}", fontsize = f_size)
ax3.plot([h0s[1]/S52g13, h0s[end]/S52g13],[1,1],color="cyan",linestyle="dashed",linewidth=4.)
# ax3.plot(h0s/S52g13, alphas[1]/gammas[13] .+ sum((Qy_gamma_alpha[:,:,2] .> Qy_gamma_alpha[:,:,3])',dims=1)'./21*((alphas[end] - alphas[1])/gammas[13]),linewidth=5.,color="orange")
ax3.pcolor(h0s/S52g13,alphas/gamma, msk', cmap="Greens",vmin=0., vmax=8.0)
ax3.text(0.1, 0.9, L"\mathbf{(b)}", fontsize=50, horizontalalignment="center",verticalalignment="center", transform = ax3.transAxes);


fig.text(-0.01, 0.5, L"S_{3/2} / S_{5/2}", ha="center", verticalalignment="center", fontsize=f_size, rotation="vertical")
fig.text(0.4275, 0.0,L"\partial_y \eta_b / S_{5/2}", verticalalignment="center", ha="center",fontsize=f_size)


savefig("./figs/CE_FCIC_slope_cont.png",bbox_inches="tight")



#########################\
PyPlot.matplotlib[:rc]("hatch",linewidth=1.5)
PyPlot.matplotlib[:rc]("hatch",color="red") # "#15B01A")

fig,ax = PyPlot.subplots(1,2,figsize=(8,3.5))
fig.tight_layout(pad=2.0)
ax2 = ax[1]; ax3 = ax[2];

f_size = 28.
l_size = 14.

ps_min = 0. #  minimum([phase_shift_12 ./ (2*pi) ; phase_shift_32 ./ (2*pi)])
ps_max = maximum([phase_shift_12 ./ (pi) ; phase_shift_32 ./ (pi)])

# fig.colorbar(pc4)
pc2=ax2.pcolormesh(h0s/S52g13,alphas/gamma,(phase_shift_12 ./ (pi))',vmin=ps_min,vmax=ps_max, cmap=PyPlot.cm.Blues,linewidth=0,rasterized=true)     # ocean_r
ax2.tick_params(labelsize=l_size)
ax2.set_ylim([alphas[1]/gamma, alphas[end]/gamma])
ax2.set_xlim([h0s[1]/S52g13, h0s[end]/S52g13])
# ax2.set_yticklabels([])
ax2.set_title(L"\Delta \phi_{3/2}", fontsize = f_size, y=1.05)
# ax2.set_xticklabels([])
ax2.plot([h0s[1]/S52g13, h0s[end]/S52g13],[1,1],color="cyan",linestyle="dashed",linewidth=2.)
ax2.contourf(h0s/S52g13,alphas/gamma, msk', hatches=["/ ", ""], alpha=0.)
# ax2.pcolor(h0s/S52g13,alphas/gamma, msk', hatch="/", alpha=0.)
ax2.text(0.1, 0.9, L"\mathbf{(a)}", fontsize=25, horizontalalignment="center",verticalalignment="center", transform = ax2.transAxes);

# ax2.plot(h0s/S52g13, alphas[1]/gammas[13] .+ sum(abs.((Qy_gamma_alpha[:,:,3] .< - Qy_gamma_alpha[:,:,2]) .- 1)',dims=1)'./21*((alphas[end] - alphas[1])/gammas[13]),linewidth=5.,color="orange")


pc3=ax3.pcolormesh(h0s/S52g13,alphas/gamma,(phase_shift_32 ./ (pi))',vmin=ps_min,vmax=ps_max, cmap=PyPlot.cm.Blues,linewidth=0,rasterized=true)
ax3.tick_params(labelsize=l_size)
ax3.set_ylim([alphas[1]/gamma, alphas[end]/gamma])
ax3.set_xlim([h0s[1]/S52g13, h0s[end]/S52g13])
ax3.set_yticklabels([])
ax3.set_title(L"\Delta \phi_{5/2}", fontsize = f_size, y=1.05)
ax3.plot([h0s[1]/S52g13, h0s[end]/S52g13],[1,1],color="cyan",linestyle="dashed",linewidth=2.)
# ax3.plot(h0s/S52g13, alphas[1]/gammas[13] .+ sum((Qy_gamma_alpha[:,:,2] .> Qy_gamma_alpha[:,:,3])',dims=1)'./21*((alphas[end] - alphas[1])/gammas[13]),linewidth=5.,color="orange")
ax3.pcolor(h0s/S52g13,alphas/gamma, msk', hatch="/", alpha=0.)
ax3.text(0.1, 0.9, L"\mathbf{(b)}", fontsize=25, horizontalalignment="center",verticalalignment="center", transform = ax3.transAxes);



fig.text(-0.02, 0.5, L"S_{3/2} / S_{5/2}", ha="center", verticalalignment="center", fontsize=f_size, rotation="vertical")
fig.text(0.45, 0.0,L"\partial_y \eta_b / S_{5/2}", verticalalignment="center", ha="center",fontsize=f_size)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cb = fig.colorbar(pc2, cax=cbar_ax)
cb.ax.tick_params(labelsize=l_size)
cb.ax.yaxis.offsetText.set(size=l_size)
cb.ax.set_yticklabels([latexstring(string.(round.(i; digits=2))) for i in cb.get_ticks()])


savefig("./figs/CE_FCIC_slope.png",bbox_inches="tight",dpi=400)

savefig("./figs/CE_FCIC_slope.pdf",bbox_inches="tight")
