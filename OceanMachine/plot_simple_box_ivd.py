import xarray as xr
import matplotlib.pyplot as plt
import cmocean
import ffmpeg

ds = xr.open_dataset("simple_box_ivd.nc")

x = ds.x_nodal.isel(z=-1) / 1000
y = ds.y_nodal.isel(z=-1) / 1000

Nt = ds.time.size

for n in range(Nt):
    u = ds.u.isel(z=-1, time=n).squeeze()
    v = ds.v.isel(z=-1, time=n).squeeze()
    η = ds.η.isel(z=-1, time=n).squeeze()
    θ = ds.θ.isel(z=-1, time=n).squeeze()

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(16, 16), dpi=200)
    plt.subplots_adjust(hspace=0.25)
    fig.suptitle(f"Simple box IVDC @ z = 0, iteration {n:02d}", fontsize=16)

    u_mesh, v_mesh = axes[0, 0], axes[0, 1]
    η_mesh, θ_mesh = axes[1, 0], axes[1, 1]

    # img_u = u.plot.pcolormesh(ax=u_mesh, vmin=-0.2, vmax=0.2, cmap=cmocean.cm.balance, add_colorbar=False)
    img_u = u_mesh.pcolormesh(x.values, y.values, u.values, vmin=-0.2, vmax=0.2, cmap=cmocean.cm.balance, shading="gouraud") 
    fig.colorbar(img_u, ax=u_mesh, label="m/s", extend="both")
    u_mesh.set_title("u-velocity")
    u_mesh.set_xlabel("x (km)")
    u_mesh.set_ylabel("z (km)")
    u_mesh.set_aspect("equal")

    # img_v = v.plot.pcolormesh(ax=v_mesh, vmin=-0.2, vmax=0.2, cmap=cmocean.cm.balance, add_colorbar=False)
    img_v = v_mesh.pcolormesh(x.values, y.values, v.values, vmin=-0.2, vmax=0.2, cmap=cmocean.cm.balance, shading="gouraud") 
    fig.colorbar(img_v, ax=v_mesh, label="m/s", extend="both")
    v_mesh.set_title("v-velocity")
    v_mesh.set_xlabel("x (km)")
    v_mesh.set_ylabel("z (km)")
    v_mesh.set_aspect("equal")

    # img_η = η.plot.pcolormesh(ax=η_mesh, vmin=-1, vmax=1, cmap=cmocean.cm.balance, add_colorbar=False)
    img_η = η_mesh.pcolormesh(x.values, y.values, η.values, vmin=-1, vmax=1, cmap=cmocean.cm.balance, shading="gouraud") 
    fig.colorbar(img_η, ax=η_mesh, label="m", extend="both")
    η_mesh.set_title("Sea surface height η")
    η_mesh.set_xlabel("x (km)")
    η_mesh.set_ylabel("z (km)")
    η_mesh.set_aspect("equal")

    # img_θ = θ.plot.pcolormesh(ax=θ_mesh, vmin=0, vmax=10, cmap=cmocean.cm.thermal, add_colorbar=False)
    img_θ = θ_mesh.pcolormesh(x.values, y.values, θ.values, vmin=0, vmax=10, cmap=cmocean.cm.thermal, shading="gouraud") 
    fig.colorbar(img_θ, ax=θ_mesh, label="°C", extend="both")
    θ_mesh.set_title("Temperature θ")
    θ_mesh.set_xlabel("x (km)")
    θ_mesh.set_ylabel("z (km)")
    θ_mesh.set_aspect("equal")

    filename = f"simple_box_ivdc_{n:02d}.png"
    print(f"Saving: {filename}")
    plt.savefig(filename)
    plt.close("all")

(
    ffmpeg
    .input(f"simple_box_ivdc_%02d.png", framerate=30)
    .output(f"simple_box_ivdc_gouraud.mp4", crf=15, pix_fmt='yuv420p')
    .overwrite_output()
    .run()
)
