<!-- Generator: GNU source-highlight 3.1.5
by Lorenzo Bettini
http://www.lorenzobettini.it
http://www.gnu.org/software/src-highlite -->
<pre><tt>

<b><font color="#0000FF">module</font></b> array_stuff
    <b><font color="#0000FF">implicit</font></b> none

<b><font color="#0000FF">contains</font></b>

    <b><font color="#0000FF">subroutine</font></b> <b><font color="#000000">read_input</font></b><font color="#990000">(</font>xyzL<font color="#990000">,</font> xyz_orig<font color="#990000">,</font> npxyz_MD<font color="#990000">,</font> ncxyz<font color="#990000">,</font> npxyz_CFD<font color="#990000">)</font>

        <font color="#009900">integer</font><font color="#990000">,</font> <b><font color="#0000FF">dimension</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">),</font> <b><font color="#0000FF">optional</font></b><font color="#990000">,</font> <b><font color="#0000FF">intent</font></b><font color="#990000">(</font>out<font color="#990000">)</font> <font color="#990000">::</font> npxyz_MD<font color="#990000">,</font> npxyz_CFD<font color="#990000">,</font> ncxyz
        <font color="#009900">double</font> <font color="#009900">precision</font><font color="#990000">,</font> <b><font color="#0000FF">dimension</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">),</font> <b><font color="#0000FF">intent</font></b><font color="#990000">(</font>out<font color="#990000">)</font>  <font color="#990000">::</font> xyzL<font color="#990000">,</font> xyz_orig

        <font color="#009900">integer</font> <font color="#990000">::</font> dummy

        <b><font color="#0000FF">open</font></b><font color="#990000">(</font><b><font color="#0000FF">unit</font></b><font color="#990000">=</font><font color="#993399">10</font><font color="#990000">,</font> <b><font color="#0000FF">file</font></b><font color="#990000">=</font><font color="#FF0000">'./fortran.in'</font><font color="#990000">)</font>
        <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> <b><font color="#000000">xyzL</font></b><font color="#990000">(</font><font color="#993399">1</font><font color="#990000">)</font>
        <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> <b><font color="#000000">xyzL</font></b><font color="#990000">(</font><font color="#993399">2</font><font color="#990000">)</font>
        <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> <b><font color="#000000">xyzL</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">)</font>
        <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> <b><font color="#000000">xyz_orig</font></b><font color="#990000">(</font><font color="#993399">1</font><font color="#990000">)</font>
        <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> <b><font color="#000000">xyz_orig</font></b><font color="#990000">(</font><font color="#993399">2</font><font color="#990000">)</font>
        <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> <b><font color="#000000">xyz_orig</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">)</font>
        <b><font color="#0000FF">if</font></b> <font color="#990000">(</font><b><font color="#000000">present</font></b><font color="#990000">(</font>npxyz_CFD<font color="#990000">))</font> <b><font color="#0000FF">then</font></b>
            <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> <b><font color="#000000">npxyz_CFD</font></b><font color="#990000">(</font><font color="#993399">1</font><font color="#990000">)</font>
            <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> <b><font color="#000000">npxyz_CFD</font></b><font color="#990000">(</font><font color="#993399">2</font><font color="#990000">)</font>
            <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> <b><font color="#000000">npxyz_CFD</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">)</font>
        <b><font color="#0000FF">else</font></b>
            <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> dummy
            <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> dummy
            <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> dummy
        <b><font color="#0000FF">endif</font></b>
        <b><font color="#0000FF">if</font></b> <font color="#990000">(</font><b><font color="#000000">present</font></b><font color="#990000">(</font>ncxyz<font color="#990000">))</font> <b><font color="#0000FF">then</font></b>
            <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> <b><font color="#000000">ncxyz</font></b><font color="#990000">(</font><font color="#993399">1</font><font color="#990000">)</font>
            <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> <b><font color="#000000">ncxyz</font></b><font color="#990000">(</font><font color="#993399">2</font><font color="#990000">)</font>
            <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> <b><font color="#000000">ncxyz</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">)</font>
        <b><font color="#0000FF">else</font></b>
            <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> dummy
            <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> dummy
            <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> dummy
        <b><font color="#0000FF">endif</font></b>
        <b><font color="#0000FF">if</font></b> <font color="#990000">(</font><b><font color="#000000">present</font></b><font color="#990000">(</font>npxyz_MD<font color="#990000">))</font> <b><font color="#0000FF">then</font></b>
            <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> <b><font color="#000000">npxyz_MD</font></b><font color="#990000">(</font><font color="#993399">1</font><font color="#990000">)</font>
            <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> <b><font color="#000000">npxyz_MD</font></b><font color="#990000">(</font><font color="#993399">2</font><font color="#990000">)</font>
            <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> <b><font color="#000000">npxyz_MD</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">)</font>
        <b><font color="#0000FF">else</font></b>
            <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> dummy
            <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> dummy
            <b><font color="#0000FF">read</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">,</font> <font color="#990000">*)</font> dummy
        <b><font color="#0000FF">endif</font></b>
        <b><font color="#0000FF">close</font></b><font color="#990000">(</font><font color="#993399">10</font><font color="#990000">)</font>

        <i><font color="#9A1900">!print*, xyzL, xyz_orig, npxyz_MD, ncxyz, npxyz_CFD</font></i>

    <b><font color="#0000FF">end</font></b> <b><font color="#0000FF">subroutine</font></b> 

    <b><font color="#0000FF">subroutine</font></b> <b><font color="#000000">fill_array</font></b><font color="#990000">(</font>A<font color="#990000">)</font>
        <b><font color="#0000FF">use</font></b> cpl<font color="#990000">,</font> <b><font color="#0000FF">only</font></b> <font color="#990000">:</font> CPL_get_olap_limits<font color="#990000">,</font> CPL_my_proc_portion<font color="#990000">,</font> <font color="#CC33CC">&amp;</font>
                        CPL_get_no_cells

        <font color="#009900">double</font> <font color="#009900">precision</font><font color="#990000">,</font> <b><font color="#0000FF">dimension</font></b><font color="#990000">(:,:,:,:),</font> <font color="#CC33CC">&amp;</font>
            <b><font color="#0000FF">allocatable</font></b><font color="#990000">,</font> <b><font color="#0000FF">intent</font></b><font color="#990000">(</font>out<font color="#990000">)</font>  <font color="#990000">::</font> A

        <b><font color="#0000FF">logical</font></b> <font color="#990000">::</font> no_error
        <font color="#009900">integer</font> <font color="#990000">::</font> i<font color="#990000">,</font>j<font color="#990000">,</font>k<font color="#990000">,</font>ii<font color="#990000">,</font>jj<font color="#990000">,</font>kk<font color="#990000">,</font>ierr
        <font color="#009900">integer</font><font color="#990000">,</font> <b><font color="#0000FF">dimension</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">)</font> <font color="#990000">::</font> Ncells
        <font color="#009900">integer</font><font color="#990000">,</font> <b><font color="#0000FF">dimension</font></b><font color="#990000">(</font><font color="#993399">6</font><font color="#990000">)</font> <font color="#990000">::</font> portion<font color="#990000">,</font> limits

        <b><font color="#0000FF">call</font></b> <b><font color="#000000">CPL_get_olap_limits</font></b><font color="#990000">(</font>limits<font color="#990000">)</font>
        <b><font color="#0000FF">call</font></b> <b><font color="#000000">CPL_my_proc_portion</font></b><font color="#990000">(</font>limits<font color="#990000">,</font> portion<font color="#990000">)</font>
        <b><font color="#0000FF">call</font></b> <b><font color="#000000">CPL_get_no_cells</font></b><font color="#990000">(</font>portion<font color="#990000">,</font> Ncells<font color="#990000">)</font>

        <b><font color="#0000FF">allocate</font></b><font color="#990000">(</font><b><font color="#000000">A</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">,</font> <b><font color="#000000">Ncells</font></b><font color="#990000">(</font><font color="#993399">1</font><font color="#990000">),</font> <b><font color="#000000">Ncells</font></b><font color="#990000">(</font><font color="#993399">2</font><font color="#990000">),</font> <b><font color="#000000">Ncells</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">)))</font>
        <b><font color="#0000FF">do</font></b> i <font color="#990000">=</font> <font color="#993399">1</font><font color="#990000">,</font><b><font color="#000000">Ncells</font></b><font color="#990000">(</font><font color="#993399">1</font><font color="#990000">)</font>
        <b><font color="#0000FF">do</font></b> j <font color="#990000">=</font> <font color="#993399">1</font><font color="#990000">,</font><b><font color="#000000">Ncells</font></b><font color="#990000">(</font><font color="#993399">2</font><font color="#990000">)</font>
        <b><font color="#0000FF">do</font></b> k <font color="#990000">=</font> <font color="#993399">1</font><font color="#990000">,</font><b><font color="#000000">Ncells</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">)</font>
            ii <font color="#990000">=</font> i <font color="#990000">+</font> <b><font color="#000000">portion</font></b><font color="#990000">(</font><font color="#993399">1</font><font color="#990000">)</font>
            jj <font color="#990000">=</font> j <font color="#990000">+</font> <b><font color="#000000">portion</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">)</font>
            kk <font color="#990000">=</font> k <font color="#990000">+</font> <b><font color="#000000">portion</font></b><font color="#990000">(</font><font color="#993399">5</font><font color="#990000">)</font>

            <b><font color="#000000">A</font></b><font color="#990000">(</font><font color="#993399">1</font><font color="#990000">,</font>i<font color="#990000">,</font>j<font color="#990000">,</font>k<font color="#990000">)</font> <font color="#990000">=</font> ii
            <b><font color="#000000">A</font></b><font color="#990000">(</font><font color="#993399">2</font><font color="#990000">,</font>i<font color="#990000">,</font>j<font color="#990000">,</font>k<font color="#990000">)</font> <font color="#990000">=</font> jj
            <b><font color="#000000">A</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">,</font>i<font color="#990000">,</font>j<font color="#990000">,</font>k<font color="#990000">)</font> <font color="#990000">=</font> kk
        <b><font color="#0000FF">enddo</font></b>
        <b><font color="#0000FF">enddo</font></b>
        <b><font color="#0000FF">enddo</font></b>

    <b><font color="#0000FF">end</font></b> <b><font color="#0000FF">subroutine</font></b> fill_array


    <b><font color="#0000FF">subroutine</font></b> <b><font color="#000000">print_array</font></b><font color="#990000">(</font>A<font color="#990000">,</font> rank<font color="#990000">,</font> no_error<font color="#990000">)</font>
        <b><font color="#0000FF">use</font></b> cpl<font color="#990000">,</font> <b><font color="#0000FF">only</font></b> <font color="#990000">:</font> CPL_get_olap_limits<font color="#990000">,</font> CPL_my_proc_portion<font color="#990000">,</font> <font color="#CC33CC">&amp;</font>
                        CPL_get_no_cells

        <b><font color="#0000FF">logical</font></b><font color="#990000">,</font> <b><font color="#0000FF">intent</font></b><font color="#990000">(</font>out<font color="#990000">),</font> <b><font color="#0000FF">optional</font></b> <font color="#990000">::</font> no_error
        <font color="#009900">integer</font><font color="#990000">,</font> <b><font color="#0000FF">intent</font></b><font color="#990000">(</font>in<font color="#990000">),</font> <b><font color="#0000FF">optional</font></b> <font color="#990000">::</font> rank
        <font color="#009900">double</font> <font color="#009900">precision</font><font color="#990000">,</font> <b><font color="#0000FF">dimension</font></b><font color="#990000">(:,:,:,:),</font> <font color="#CC33CC">&amp;</font>
            <b><font color="#0000FF">allocatable</font></b><font color="#990000">,</font> <b><font color="#0000FF">intent</font></b><font color="#990000">(</font>in<font color="#990000">)</font>  <font color="#990000">::</font> A

        <font color="#009900">integer</font> <font color="#990000">::</font> i<font color="#990000">,</font>j<font color="#990000">,</font>k<font color="#990000">,</font>ii<font color="#990000">,</font>jj<font color="#990000">,</font>kk<font color="#990000">,</font>ierr
        <font color="#009900">integer</font><font color="#990000">,</font> <b><font color="#0000FF">dimension</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">)</font> <font color="#990000">::</font> Ncells
        <font color="#009900">integer</font><font color="#990000">,</font> <b><font color="#0000FF">dimension</font></b><font color="#990000">(</font><font color="#993399">6</font><font color="#990000">)</font> <font color="#990000">::</font> portion<font color="#990000">,</font> limits

        <b><font color="#0000FF">call</font></b> <b><font color="#000000">CPL_get_olap_limits</font></b><font color="#990000">(</font>limits<font color="#990000">)</font>
        <b><font color="#0000FF">call</font></b> <b><font color="#000000">CPL_my_proc_portion</font></b><font color="#990000">(</font>limits<font color="#990000">,</font> portion<font color="#990000">)</font>
        <b><font color="#0000FF">call</font></b> <b><font color="#000000">CPL_get_no_cells</font></b><font color="#990000">(</font>portion<font color="#990000">,</font> Ncells<font color="#990000">)</font>

        <b><font color="#0000FF">if</font></b> <font color="#990000">(</font><b><font color="#000000">any</font></b><font color="#990000">(</font>portion <font color="#990000">.</font>ne<font color="#990000">.</font> <font color="#993399">-666</font><font color="#990000">))</font> <b><font color="#0000FF">then</font></b>
            no_error <font color="#990000">=</font> <font color="#990000">.</font>true<font color="#990000">.</font>
            <b><font color="#0000FF">do</font></b> i <font color="#990000">=</font> <font color="#993399">1</font><font color="#990000">,</font><b><font color="#000000">Ncells</font></b><font color="#990000">(</font><font color="#993399">1</font><font color="#990000">)</font>
            <b><font color="#0000FF">do</font></b> j <font color="#990000">=</font> <font color="#993399">1</font><font color="#990000">,</font><b><font color="#000000">Ncells</font></b><font color="#990000">(</font><font color="#993399">2</font><font color="#990000">)</font>
            <b><font color="#0000FF">do</font></b> k <font color="#990000">=</font> <font color="#993399">1</font><font color="#990000">,</font><b><font color="#000000">Ncells</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">)</font>
                ii <font color="#990000">=</font> i <font color="#990000">+</font> <b><font color="#000000">portion</font></b><font color="#990000">(</font><font color="#993399">1</font><font color="#990000">)</font>
                jj <font color="#990000">=</font> j <font color="#990000">+</font> <b><font color="#000000">portion</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">)</font>
                kk <font color="#990000">=</font> k <font color="#990000">+</font> <b><font color="#000000">portion</font></b><font color="#990000">(</font><font color="#993399">5</font><font color="#990000">)</font>

                <b><font color="#0000FF">if</font></b> <font color="#990000">((</font><b><font color="#000000">dble</font></b><font color="#990000">(</font>ii<font color="#990000">)</font> <font color="#990000">-</font> <b><font color="#000000">A</font></b><font color="#990000">(</font><font color="#993399">1</font><font color="#990000">,</font>i<font color="#990000">,</font>j<font color="#990000">,</font>k<font color="#990000">))</font> <font color="#990000">.</font>gt<font color="#990000">.</font> <font color="#993399">1e-8</font><font color="#990000">)</font> <b><font color="#0000FF">then</font></b> 
                    <b><font color="#0000FF">print</font></b><font color="#FF0000">'(a,2i5,a,i5,a,i6,a,f10.5)'</font><font color="#990000">,</font> <font color="#FF0000">"ERROR -- portion: "</font><font color="#990000">,</font> <b><font color="#000000">portion</font></b><font color="#990000">(</font><font color="#993399">1</font><font color="#990000">:</font><font color="#993399">2</font><font color="#990000">),</font> <font color="#CC33CC">&amp;</font> 
                           <font color="#FF0000">" MD rank: "</font><font color="#990000">,</font> rank<font color="#990000">,</font> <font color="#FF0000">" cell i: "</font><font color="#990000">,</font>ii<font color="#990000">,</font> <font color="#CC33CC">&amp;</font> 
                           <font color="#FF0000">" recv_array: "</font><font color="#990000">,</font> <b><font color="#000000">A</font></b><font color="#990000">(</font><font color="#993399">1</font><font color="#990000">,</font>i<font color="#990000">,</font>j<font color="#990000">,</font>k<font color="#990000">)</font>
                    no_error <font color="#990000">=</font> <font color="#990000">.</font>false<font color="#990000">.</font>
                <b><font color="#0000FF">endif</font></b>
                <b><font color="#0000FF">if</font></b> <font color="#990000">((</font><b><font color="#000000">dble</font></b><font color="#990000">(</font>jj<font color="#990000">)</font> <font color="#990000">-</font> <b><font color="#000000">A</font></b><font color="#990000">(</font><font color="#993399">2</font><font color="#990000">,</font>i<font color="#990000">,</font>j<font color="#990000">,</font>k<font color="#990000">))</font> <font color="#990000">.</font>gt<font color="#990000">.</font> <font color="#993399">1e-8</font><font color="#990000">)</font> <b><font color="#0000FF">then</font></b> 
                    <b><font color="#0000FF">print</font></b><font color="#FF0000">'(a,2i5,a,i5,a,i6,a,f10.5)'</font><font color="#990000">,</font> <font color="#FF0000">"ERROR -- portion: "</font><font color="#990000">,</font> <b><font color="#000000">portion</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">:</font><font color="#993399">4</font><font color="#990000">),</font> <font color="#CC33CC">&amp;</font> 
                           <font color="#FF0000">" MD rank: "</font><font color="#990000">,</font> rank<font color="#990000">,</font> <font color="#FF0000">" cell j: "</font><font color="#990000">,</font> jj <font color="#990000">,</font> <font color="#CC33CC">&amp;</font> 
                           <font color="#FF0000">" recv_array: "</font><font color="#990000">,</font> <b><font color="#000000">A</font></b><font color="#990000">(</font><font color="#993399">2</font><font color="#990000">,</font>i<font color="#990000">,</font>j<font color="#990000">,</font>k<font color="#990000">)</font>
                    no_error <font color="#990000">=</font> <font color="#990000">.</font>false<font color="#990000">.</font>  
                <b><font color="#0000FF">endif</font></b>
                <b><font color="#0000FF">if</font></b> <font color="#990000">((</font><b><font color="#000000">dble</font></b><font color="#990000">(</font>kk<font color="#990000">)</font> <font color="#990000">-</font> <b><font color="#000000">A</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">,</font>i<font color="#990000">,</font>j<font color="#990000">,</font>k<font color="#990000">))</font> <font color="#990000">.</font>gt<font color="#990000">.</font> <font color="#993399">1e-8</font><font color="#990000">)</font> <b><font color="#0000FF">then</font></b> 
                    <b><font color="#0000FF">print</font></b><font color="#FF0000">'(a,2i5,a,i5,a,i6,a,f10.5)'</font><font color="#990000">,</font> <font color="#FF0000">"ERROR -- portion: "</font><font color="#990000">,</font> <b><font color="#000000">portion</font></b><font color="#990000">(</font><font color="#993399">5</font><font color="#990000">:</font><font color="#993399">6</font><font color="#990000">),</font> <font color="#CC33CC">&amp;</font> 
                           <font color="#FF0000">" MD rank: "</font><font color="#990000">,</font> rank<font color="#990000">,</font> <font color="#FF0000">" cell k: "</font><font color="#990000">,</font> kk <font color="#990000">,</font> <font color="#CC33CC">&amp;</font> 
                           <font color="#FF0000">" recv_array: "</font><font color="#990000">,</font> <b><font color="#000000">A</font></b><font color="#990000">(</font><font color="#993399">3</font><font color="#990000">,</font>i<font color="#990000">,</font>j<font color="#990000">,</font>k<font color="#990000">)</font>
                    no_error <font color="#990000">=</font> <font color="#990000">.</font>false<font color="#990000">.</font>
                <b><font color="#0000FF">endif</font></b>
            <b><font color="#0000FF">enddo</font></b>
            <b><font color="#0000FF">enddo</font></b>
            <b><font color="#0000FF">enddo</font></b>

        <b><font color="#0000FF">endif</font></b>

        <b><font color="#0000FF">if</font></b> <font color="#990000">(</font>no_error<font color="#990000">)</font> <b><font color="#0000FF">print</font></b><font color="#990000">*,</font> <font color="#FF0000">"CFD -- NO PROBLEMS HAVE BEEN DETECTED!"</font>

    <b><font color="#0000FF">end</font></b> <b><font color="#0000FF">subroutine</font></b> print_array


<b><font color="#0000FF">end</font></b> <b><font color="#0000FF">module</font></b> array_stuff

</tt></pre>
