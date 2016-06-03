<!-- Generator: GNU source-highlight 3.1.5
by Lorenzo Bettini
http://www.lorenzobettini.it
http://www.gnu.org/software/src-highlite -->
<pre><tt>

<b><font color="#0000FF">program</font></b> cfd_cpl_example
    <b><font color="#0000FF">use</font></b> cpl<font color="#990000">,</font> <b><font color="#0000FF">only</font></b> <font color="#990000">:</font> CPL_init<font color="#990000">,</font> CPL_setup_cfd<font color="#990000">,</font> <font color="#CC33CC">&amp;</font> 
                    CPL_get_olap_limits<font color="#990000">,</font> CPL_my_proc_portion<font color="#990000">,</font> <font color="#CC33CC">&amp;</font>
                    CPL_get_no_cells<font color="#990000">,</font> CPL_send<font color="#990000">,</font> CPL_recv
    <b><font color="#0000FF">use</font></b> array_stuff<font color="#990000">,</font> <b><font color="#0000FF">only</font></b> <font color="#990000">:</font> fill_array<font color="#990000">,</font> print_array
    <b><font color="#0000FF">use</font></b> mpi
    <b><font color="#0000FF">implicit</font></b> none

    <b><font color="#0000FF">logical</font></b> <font color="#990000">::</font> recv_flag<font color="#990000">,</font> send_flag<font color="#990000">,</font> NO_ERROR
    <font color="#009900">integer</font> <font color="#990000">::</font> i<font color="#990000">,</font>j<font color="#990000">,</font>k<font color="#990000">,</font>ii<font color="#990000">,</font>jj<font color="#990000">,</font>kk<font color="#990000">,</font>ierr<font color="#990000">,</font>errorcode
    <font color="#009900">integer</font> <font color="#990000">::</font> NPx<font color="#990000">,</font> NPy<font color="#990000">,</font> NPz<font color="#990000">,</font> NProcs<font color="#990000">,</font> rank
    <font color="#009900">integer</font> <font color="#990000">::</font> nprocs_realm<font color="#990000">,</font> nsteps<font color="#990000">,</font> initialstep
    <font color="#009900">integer</font> <font color="#990000">::</font> CART_COMM<font color="#990000">,</font> CFD_COMM
    <font color="#009900">integer</font><font color="#990000">,</font> <b><font color="#0000FF">parameter</font></b> <font color="#990000">::</font> cfd_realm<font color="#990000">=</font><font color="#993399">1</font>
    <font color="#009900">integer</font><font color="#990000">,</font> <b><font color="#0000FF">dimension</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">)</font> <font color="#990000">::</font> npxyz<font color="#990000">,</font> Ncells<font color="#990000">,</font> ncxyz
    <font color="#009900">integer</font><font color="#990000">,</font> <b><font color="#0000FF">dimension</font></b><font color="#990000">(</font><font color="#993399">6</font><font color="#990000">)</font> <font color="#990000">::</font> portion<font color="#990000">,</font> limits
    <font color="#009900">double</font> <font color="#009900">precision</font> <font color="#990000">::</font> dt<font color="#990000">,</font> density
    <font color="#009900">double</font> <font color="#009900">precision</font><font color="#990000">,</font> <b><font color="#0000FF">dimension</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">)</font>  <font color="#990000">::</font> xyzL<font color="#990000">,</font> xyz_orig
    <font color="#009900">double</font> <font color="#009900">precision</font><font color="#990000">,</font> <b><font color="#0000FF">dimension</font></b><font color="#990000">(:,:,:,:),</font> <font color="#CC33CC">&amp;</font> 
         <b><font color="#0000FF">allocatable</font></b>  <font color="#990000">::</font> send_array<font color="#990000">,</font> recv_array

    <i><font color="#9A1900">!Initialise MPI</font></i>
    <b><font color="#0000FF">call</font></b> <b><font color="#000000">MPI_Init</font></b><font color="#990000">(</font>ierr<font color="#990000">)</font>

    <i><font color="#9A1900">!Create MD Comm by spliting world</font></i>
    <b><font color="#0000FF">call</font></b> <b><font color="#000000">CPL_init</font></b><font color="#990000">(</font>CFD_realm<font color="#990000">,</font> CFD_COMM<font color="#990000">,</font> ierr<font color="#990000">)</font>

    <i><font color="#9A1900">!Parameters</font></i>
    dt <font color="#990000">=</font> <font color="#993399">0.1</font>
    density <font color="#990000">=</font> <font color="#993399">0.8</font>
    nsteps <font color="#990000">=</font> <font color="#993399">100</font>

    <i><font color="#9A1900">! Parameters of the cpu topology (cartesian grid)</font></i>
    xyzL <font color="#990000">=</font> <font color="#990000">(/</font><font color="#993399">10.d0</font><font color="#990000">,</font> <font color="#993399">10.d0</font><font color="#990000">,</font> <font color="#993399">10.d0</font><font color="#990000">/)</font>
    xyz_orig <font color="#990000">=</font> <font color="#990000">(/</font><font color="#993399">0.d0</font><font color="#990000">,</font> <font color="#993399">0.d0</font><font color="#990000">,</font> <font color="#993399">0.d0</font><font color="#990000">/)</font>
    npxyz <font color="#990000">=</font> <font color="#990000">(/</font> <font color="#993399">2</font><font color="#990000">,</font> <font color="#993399">2</font><font color="#990000">,</font> <font color="#993399">1</font><font color="#990000">/)</font>
    ncxyz <font color="#990000">=</font> <font color="#990000">(/</font> <font color="#993399">64</font><font color="#990000">,</font> <font color="#993399">18</font><font color="#990000">,</font> <font color="#993399">64</font> <font color="#990000">/)</font>

    <i><font color="#9A1900">!call read_input(xyzL=xyzL, xyz_orig=xyz_orig, &amp; </font></i>
    <i><font color="#9A1900">!                npxyz_CFD=npxyz, ncxyz=ncxyz)</font></i>

    <i><font color="#9A1900">! Create communicators and check that number of processors is consistent</font></i>
	<b><font color="#0000FF">call</font></b> <b><font color="#000000">MPI_Comm_size</font></b><font color="#990000">(</font>CFD_COMM<font color="#990000">,</font> nprocs_realm<font color="#990000">,</font> ierr<font color="#990000">)</font> 
    <b><font color="#0000FF">if</font></b> <font color="#990000">(</font>nprocs_realm <font color="#990000">.</font>ne<font color="#990000">.</font> <b><font color="#000000">product</font></b><font color="#990000">(</font>npxyz<font color="#990000">))</font> <b><font color="#0000FF">then</font></b>
        <b><font color="#0000FF">print</font></b><font color="#990000">*,</font> <font color="#FF0000">"Non-coherent number of processes"</font>
        <b><font color="#0000FF">call</font></b> <b><font color="#000000">MPI_Abort</font></b><font color="#990000">(</font>MPI_COMM_WORLD<font color="#990000">,</font> errorcode<font color="#990000">,</font> ierr<font color="#990000">)</font>
    <b><font color="#0000FF">endif</font></b>

    <i><font color="#9A1900">!Setup cartesian topology</font></i>
    <b><font color="#0000FF">call</font></b> <b><font color="#000000">MPI_comm_rank</font></b><font color="#990000">(</font>CFD_COMM<font color="#990000">,</font> rank<font color="#990000">,</font> ierr<font color="#990000">)</font>
	<b><font color="#0000FF">call</font></b> <b><font color="#000000">MPI_Cart_create</font></b><font color="#990000">(</font>CFD_COMM<font color="#990000">,</font> <font color="#993399">3</font><font color="#990000">,</font> npxyz<font color="#990000">,</font> <font color="#990000">(/.</font>true<font color="#990000">.,.</font>true<font color="#990000">.,.</font>true<font color="#990000">./),</font> <font color="#CC33CC">&amp;</font> 
                         <font color="#990000">.</font>true<font color="#990000">.,</font> CART_COMM<font color="#990000">,</font> ierr<font color="#990000">)</font>

    <i><font color="#9A1900">!Coupler setup</font></i>
    <b><font color="#0000FF">call</font></b> <b><font color="#000000">CPL_setup_cfd</font></b><font color="#990000">(</font>nsteps<font color="#990000">,</font> dt<font color="#990000">,</font> CART_COMM<font color="#990000">,</font> xyzL<font color="#990000">,</font> xyz_orig<font color="#990000">,</font> ncxyz<font color="#990000">,</font> density<font color="#990000">)</font>

    <i><font color="#9A1900">!Get detail for grid</font></i>
    <b><font color="#0000FF">call</font></b> <b><font color="#000000">CPL_get_olap_limits</font></b><font color="#990000">(</font>limits<font color="#990000">)</font>
    <b><font color="#0000FF">call</font></b> <b><font color="#000000">CPL_my_proc_portion</font></b><font color="#990000">(</font>limits<font color="#990000">,</font> portion<font color="#990000">)</font>
    <b><font color="#0000FF">call</font></b> <b><font color="#000000">CPL_get_no_cells</font></b><font color="#990000">(</font>portion<font color="#990000">,</font> Ncells<font color="#990000">)</font>

    <i><font color="#9A1900">!Pack array with cell data</font></i>
    <b><font color="#0000FF">call</font></b> <b><font color="#000000">fill_array</font></b><font color="#990000">(</font>send_array<font color="#990000">)</font>
    <b><font color="#0000FF">call</font></b> <b><font color="#000000">CPL_send</font></b><font color="#990000">(</font>send_array<font color="#990000">,</font> limits<font color="#990000">,</font> send_flag<font color="#990000">)</font>

    <i><font color="#9A1900">!Block before checking if successful</font></i>
    <b><font color="#0000FF">call</font></b> <b><font color="#000000">MPI_Barrier</font></b><font color="#990000">(</font>MPI_COMM_WORLD<font color="#990000">,</font> ierr<font color="#990000">)</font>

    <i><font color="#9A1900">!Coupled Recieve and print</font></i>
    <b><font color="#0000FF">allocate</font></b><font color="#990000">(</font><b><font color="#000000">recv_array</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">,</font> <b><font color="#000000">Ncells</font></b><font color="#990000">(</font><font color="#993399">1</font><font color="#990000">),</font> <b><font color="#000000">Ncells</font></b><font color="#990000">(</font><font color="#993399">2</font><font color="#990000">),</font> <b><font color="#000000">Ncells</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">)))</font>
    recv_array <font color="#990000">=</font> <font color="#993399">0.d0</font>
    <b><font color="#0000FF">call</font></b> <b><font color="#000000">CPL_recv</font></b><font color="#990000">(</font>recv_array<font color="#990000">,</font> limits<font color="#990000">,</font> recv_flag<font color="#990000">)</font>
    <b><font color="#0000FF">call</font></b> <b><font color="#000000">print_array</font></b><font color="#990000">(</font>recv_array<font color="#990000">,</font> rank<font color="#990000">,</font> no_error<font color="#990000">)</font>

    <i><font color="#9A1900">!Block before checking if successful</font></i>
    <b><font color="#0000FF">call</font></b> <b><font color="#000000">MPI_Barrier</font></b><font color="#990000">(</font>MPI_COMM_WORLD<font color="#990000">,</font> ierr<font color="#990000">)</font>
    <b><font color="#0000FF">if</font></b> <font color="#990000">(</font>no_error<font color="#990000">)</font> <b><font color="#0000FF">print</font></b><font color="#990000">*,</font> <font color="#FF0000">"CFD -- NO PROBLEMS HAVE BEEN DETECTED!"</font>

<b><font color="#0000FF">end</font></b> <b><font color="#0000FF">program</font></b> cfd_cpl_example
</tt></pre>
