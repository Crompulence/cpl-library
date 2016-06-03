<!-- Generator: GNU source-highlight 3.1.5
by Lorenzo Bettini
http://www.lorenzobettini.it
http://www.gnu.org/software/src-highlite -->
<pre><tt><b><font color="#000080">from</font></b> mpi4py <b><font color="#000080">import</font></b> MPI
<b><font color="#000080">from</font></b> cplpy<font color="#990000">.</font>cpl <b><font color="#000080">import</font></b> CPL 
<b><font color="#000080">import</font></b> numpy as np

CPL <font color="#990000">=</font> <b><font color="#000000">CPL</font></b><font color="#990000">()</font>
nsteps <font color="#990000">=</font> <font color="#993399">1</font>
dt <font color="#990000">=</font> <font color="#993399">0.2</font>
density <font color="#990000">=</font> <font color="#993399">0.8</font>

<i><font color="#9A1900"># Parameters of the cpu topology (cartesian grid)</font></i>
NPx <font color="#990000">=</font> <font color="#993399">2</font>
NPy <font color="#990000">=</font> <font color="#993399">2</font>
NPz <font color="#990000">=</font> <font color="#993399">1</font>
NProcs <font color="#990000">=</font> NPx<font color="#990000">*</font>NPy<font color="#990000">*</font>NPz

<i><font color="#9A1900"># Parameters of the mesh topology (cartesian grid)</font></i>
ncxyz <font color="#990000">=</font> np<font color="#990000">.</font><b><font color="#000000">array</font></b><font color="#990000">([</font><font color="#993399">64</font><font color="#990000">,</font> <font color="#993399">18</font><font color="#990000">,</font> <font color="#993399">64</font><font color="#990000">],</font> order<font color="#990000">=</font><font color="#FF0000">'F'</font><font color="#990000">,</font> dtype<font color="#990000">=</font>np<font color="#990000">.</font>int32<font color="#990000">)</font>
xyzL <font color="#990000">=</font> np<font color="#990000">.</font><b><font color="#000000">array</font></b><font color="#990000">([</font><font color="#993399">10.0</font><font color="#990000">,</font> <font color="#993399">10.0</font><font color="#990000">,</font> <font color="#993399">10.0</font><font color="#990000">],</font> order<font color="#990000">=</font><font color="#FF0000">'F'</font><font color="#990000">,</font> dtype<font color="#990000">=</font>np<font color="#990000">.</font>float64<font color="#990000">)</font>
xyz_orig <font color="#990000">=</font> np<font color="#990000">.</font><b><font color="#000000">array</font></b><font color="#990000">([</font><font color="#993399">0.0</font><font color="#990000">,</font> <font color="#993399">0.0</font><font color="#990000">,</font> <font color="#993399">0.0</font><font color="#990000">],</font> order<font color="#990000">=</font><font color="#FF0000">'F'</font><font color="#990000">,</font> dtype<font color="#990000">=</font>np<font color="#990000">.</font>float64<font color="#990000">)</font>

<i><font color="#9A1900"># Create communicators and check that number of processors is consistent</font></i>
realm_comm <font color="#990000">=</font> CPL<font color="#990000">.</font><b><font color="#000000">init</font></b><font color="#990000">(</font>CPL<font color="#990000">.</font>CFD_REALM<font color="#990000">)</font>
nprocs_realm <font color="#990000">=</font> realm_comm<font color="#990000">.</font><b><font color="#000000">Get_size</font></b><font color="#990000">()</font>

<b><font color="#0000FF">if</font></b> <font color="#990000">(</font>nprocs_realm <font color="#990000">!=</font> NProcs<font color="#990000">):</font>
    <b><font color="#0000FF">print</font></b> <font color="#FF0000">"ERROR: Non-coherent number of processors."</font>
    MPI<font color="#990000">.</font><b><font color="#000000">Abort</font></b><font color="#990000">(</font>errorcode<font color="#990000">=</font><font color="#993399">1</font><font color="#990000">)</font>

cart_comm <font color="#990000">=</font> realm_comm<font color="#990000">.</font><b><font color="#000000">Create_cart</font></b><font color="#990000">([</font>NPx<font color="#990000">,</font> NPy<font color="#990000">,</font> NPz<font color="#990000">])</font>
CPL<font color="#990000">.</font><b><font color="#000000">setup_cfd</font></b><font color="#990000">(</font>nsteps<font color="#990000">,</font> dt<font color="#990000">,</font> cart_comm<font color="#990000">,</font> xyzL<font color="#990000">,</font> xyz_orig<font color="#990000">,</font> ncxyz<font color="#990000">,</font> density<font color="#990000">)</font>

cart_rank <font color="#990000">=</font> cart_comm<font color="#990000">.</font><b><font color="#000000">Get_rank</font></b><font color="#990000">()</font>
limits <font color="#990000">=</font> CPL<font color="#990000">.</font><b><font color="#000000">get_olap_limits</font></b><font color="#990000">()</font>
portion <font color="#990000">=</font> CPL<font color="#990000">.</font><b><font color="#000000">my_proc_portion</font></b><font color="#990000">(</font>limits<font color="#990000">)</font>
<font color="#990000">[</font>ncxl<font color="#990000">,</font> ncyl<font color="#990000">,</font> nczl<font color="#990000">]</font> <font color="#990000">=</font> CPL<font color="#990000">.</font><b><font color="#000000">get_no_cells</font></b><font color="#990000">(</font>portion<font color="#990000">)</font>
send_array <font color="#990000">=</font> np<font color="#990000">.</font><b><font color="#000000">ones</font></b><font color="#990000">((</font><font color="#993399">3</font><font color="#990000">,</font> ncxl<font color="#990000">,</font> ncyl<font color="#990000">,</font> nczl<font color="#990000">),</font> order<font color="#990000">=</font><font color="#FF0000">'F'</font><font color="#990000">,</font> dtype<font color="#990000">=</font>np<font color="#990000">.</font>float64<font color="#990000">)</font>

<b><font color="#0000FF">for</font></b> i <b><font color="#0000FF">in</font></b> <b><font color="#000000">range</font></b><font color="#990000">(</font>portion<font color="#990000">[</font><font color="#993399">0</font><font color="#990000">],</font>portion<font color="#990000">[</font><font color="#993399">1</font><font color="#990000">]+</font><font color="#993399">1</font><font color="#990000">):</font>
    <b><font color="#0000FF">for</font></b> j <b><font color="#0000FF">in</font></b> <b><font color="#000000">range</font></b><font color="#990000">(</font>portion<font color="#990000">[</font><font color="#993399">2</font><font color="#990000">],</font>portion<font color="#990000">[</font><font color="#993399">3</font><font color="#990000">]+</font><font color="#993399">1</font><font color="#990000">):</font>
        <b><font color="#0000FF">for</font></b> k <b><font color="#0000FF">in</font></b> <b><font color="#000000">range</font></b><font color="#990000">(</font>portion<font color="#990000">[</font><font color="#993399">4</font><font color="#990000">],</font>portion<font color="#990000">[</font><font color="#993399">5</font><font color="#990000">]+</font><font color="#993399">1</font><font color="#990000">):</font>
            ii <font color="#990000">=</font> i <font color="#990000">-</font> portion<font color="#990000">[</font><font color="#993399">0</font><font color="#990000">]</font>
            jj <font color="#990000">=</font> j <font color="#990000">-</font> portion<font color="#990000">[</font><font color="#993399">2</font><font color="#990000">]</font>
            kk <font color="#990000">=</font> k <font color="#990000">-</font> portion<font color="#990000">[</font><font color="#993399">4</font><font color="#990000">]</font>

            send_array<font color="#990000">[</font><font color="#993399">0</font><font color="#990000">,</font>ii<font color="#990000">,</font>jj<font color="#990000">,</font>kk<font color="#990000">]</font> <font color="#990000">=</font> i<font color="#990000">+</font><font color="#993399">1</font>
            send_array<font color="#990000">[</font><font color="#993399">1</font><font color="#990000">,</font>ii<font color="#990000">,</font>jj<font color="#990000">,</font>kk<font color="#990000">]</font> <font color="#990000">=</font> j<font color="#990000">+</font><font color="#993399">1</font>
            send_array<font color="#990000">[</font><font color="#993399">2</font><font color="#990000">,</font>ii<font color="#990000">,</font>jj<font color="#990000">,</font>kk<font color="#990000">]</font> <font color="#990000">=</font> k<font color="#990000">+</font><font color="#993399">1</font>

ierr <font color="#990000">=</font> CPL<font color="#990000">.</font><b><font color="#000000">send</font></b><font color="#990000">(</font>send_array<font color="#990000">,</font> limits<font color="#990000">)</font>

MPI<font color="#990000">.</font>COMM_WORLD<font color="#990000">.</font><b><font color="#000000">Barrier</font></b><font color="#990000">()</font>
</tt></pre>
