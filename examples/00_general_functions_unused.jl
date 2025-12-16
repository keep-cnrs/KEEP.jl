using Logging

function should_verbose()
    return Logging.min_enabled_level(current_logger()) <= Info
end