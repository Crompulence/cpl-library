! Top level module containing every the setup and main coupler routines
module CPL
    use coupler_module
    use coupler
    use coupler_write
#ifdef JSON_SUPPORT
    use io
#endif 
end module CPL
