function CLIMA_slicer(solver_config)
    ## Slicing routine for structured CLIMA meshes
    # Define number nodal points along each axis of an element and number of elements along each
    # axis of topology
    Np=[4,4,4]
    Nn=[Np[1]+1,Np[2]+1,Np[3]+1]
    Ne=[20,20,20]

    # Get X, Y, Z coords and make sure they are in a CPU array to make life easier.
    Gnd=reshape(solver_config.dg.grid.vgeo,(Nn[1],Nn[2],Nn[3],16,Ne[1],Ne[2],Ne[3]))
    xvals_raw=Gnd[:,:,:,13,:,:,:]
    yvals_raw=Gnd[:,:,:,14,:,:,:]
    zvals_raw=Gnd[:,:,:,15,:,:,:]
    xvals_cpu=ones(size(xvals_raw))
    yvals_cpu=ones(size(yvals_raw))
    zvals_cpu=ones(size(zvals_raw))
    copyto!(xvals_cpu,xvals_raw)
    copyto!(yvals_cpu,yvals_raw)
    copyto!(zvals_cpu,zvals_raw)

    # Get field of interest and put on CPU too
    Qnd=reshape(solver_config.Q.realdata,(Nn[1],Nn[2],Nn[3],4,Ne[1],Ne[2],Ne[3]));
    iFld=1 # u      ( m/s)
    iFld=2 # v      ( m/s)
    iFld=3 # Eta    (   m)
    iFld=4 # Theta  (  oC)
    phi1_raw=Qnd[:,:,:,iFld,:,:,:]
    phi1_cpu=ones(size(phi1_raw))
    copyto!(phi1_cpu,phi1_raw)
 
    # Have coords and values in CPU arrays
    # Now stick things in a DataFrame to make manipulation easier
    df = DataFrame(X=xvals_cpu[:],
                   Y=yvals_cpu[:],
                   Z=zvals_cpu[:],
                   F=phi1_cpu[:])
    println(" ")
    println( first(df,6) )
    println(" ")

    # Get underlying coords by finding unique vals to some tolerance
    # Check lengths ( should be Np*Ne+1 )
    deps=1.e-2
    dedup(x,d)=(
     # x must be sorted before this
     for i=2:length(x); dx=x[i]-x[i-1];
      if abs(dx) < d;
       x[i]=x[i-1];
      end;
     end;
     return unique(x);
    )
    lx_exp=Np[1]*Ne[1]+1
    xlocs=sort!( unique(df,:X), [:X]              )[!,:X]
    xlocs=dedup(xlocs,deps)
    println(" ")
    println("xlocs = ",xlocs)
    println(" ")
    @assert length(xlocs) == lx_exp
    ly_exp=Np[2]*Ne[2]+1
    ylocs=sort!( unique(df,:Y), [:Y]              )[!,:Y]
    ylocs=dedup(ylocs,deps)
    println(" ")
    println("ylocs = ",ylocs)
    println(" ")
    @assert length(ylocs) == ly_exp
    lz_exp=Np[3]*Ne[3]+1
    zlocs=sort!( unique(df,:Z), [:Z], rev=(false) )[!,:Z]
    zlocs=dedup(zlocs,deps)
    println(" ")
    println("zlocs = ",zlocs)
    println(" ")
    @assert length(zlocs) == lz_exp

    # Now get slices by finding points that fit mesh within tolerance geps
    geps=1.e-5

    # Example 1. 
    # Line in Z with Z descending retaining jump points.
    line_z_val=df[ ( isapprox.( df.X , xlocs[1];atol=geps) ) .& 
		   ( isapprox.( df.Y , ylocs[1];atol=geps) ), : ]
    println(" ")
    println(" line_z_val")
    println( sort!(line_z_val, [:Z], rev=(true) ) )
    println(" ")
   
    s_choice="YZ"

    # XZ slice with X fastest and plot 2d heatmap
    # Set slicing params
    if s_choice == "XZ"
     const_ax  = df.Y
     const_val = ylocs[5]
     xc=xlocs
     yc=zlocs
     slice_2d_val=df[ ( isapprox.( const_ax , const_val;atol=geps) ), : ]
     ax1       = slice_2d_val.X
     xsym      = :X
     ax2       = slice_2d_val.Z
     ysym      = :Z
    end

    # YZ slice with Y fastest and plot 2d heatmap
    # Set slicing params
    if s_choice == "YZ"
     const_ax  = df.X
     const_val = xlocs[5]
     xc=ylocs
     yc=zlocs
     slice_2d_val=df[ ( isapprox.( const_ax , const_val;atol=geps) ), : ]
     ax1       = slice_2d_val.Y
     xsym      = :Y
     ax2       = slice_2d_val.Z
     ysym      = :Z
    end
     
    # XY slice with X fastest and plot 2d heatmap
    # Set slicing params
    if s_choice == "XY"
     const_ax  = df.Z
     const_val = zlocs[end]
     xc=xlocs
     yc=ylocs
     slice_2d_val=df[ ( isapprox.( const_ax , const_val;atol=geps) ), : ]
     ax1       = slice_2d_val.X
     xsym      = :X
     ax2       = slice_2d_val.Y
     ysym      = :Y
    end

    # Now have a slice specified and with points in native order
    # whatever that may be
    
    # Normalize column and row coords to fit standard locs in xc and yc
    # - this helps in dealing with jump points, and is also a check 
    # - that things are working.
    # x-dimension of slice
    i=0
    for loc in ax1
	    l_diff=loc.-xc
	    i=i+1
	    v,j=findmin(abs.(l_diff))
	    # if we understood grid coorectly all points should be within
	    # geps of a standard loc coordinate value.
	    @assert v <= geps
	    if v != 0
              slice_2d_val[i,xsym]=xc[j]
            end
    end
    # y-dimension of slice
    i=0
    for loc in ax2
	    l_diff=loc.-yc
      	    i=i+1
	    v,j=findmin(abs.(l_diff))
	    # if we understood grid coorectly all points should be within
	    # geps of a standard loc coordinate value.
	    @assert v <= geps
	    if v != 0
              slice_2d_val[i,ysym]=yc[j]
            end
    end

    println(" ")
    println(" slice_2d_val")
    println( first( sort(slice_2d_val, [xsym,ysym], rev=(false,false) ), 400) )
    println(" ")

      # Make a scatter plot - include "jump" points
      # TBD

      # Make a contourf plot - filtering "jump" points.
      # Here we filter by selecting one of them.
      # xc=xlocs
      lx=length(xc)
      # yc=zlocs
      ly=length(yc)
      # Take first jump point at element edges
      # ( Needed for a structured plot since xc and yc are Np*Ne+1, whereas )
      # ( slice_xz_val will be (Np+1)*Ne                                    )
      slice_2d_ct=unique(slice_2d_val,[xsym,ysym])
      # Now need to sort so valuees are ordered with coord
      sort!(slice_2d_ct,[xsym,ysym])
      println( size(slice_2d_ct) )
      zvals=reshape(slice_2d_ct[!,:F],(lx,ly))
      println( lx )
      println( ly )
      # savefig( contour(yc,xc,zvals), "fooo.png" )
      savefig( heatmap(xc,yc,zvals,c=:gist_ncar), "fooo.png" )

      # Save Dataframe to csv - include "jump" points.

   exit()
     
end

