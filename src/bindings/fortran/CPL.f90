! Top level module containing every the setup and main coupler routines
module CPL
    use coupler_module
    use coupler
#ifdef JSON_SUPPORT
    use io
#endif 
end module CPL
