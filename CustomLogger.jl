using Logging

debuglogger = ConsoleLogger(stderr, Logging.Debug)

global_logger(debuglogger)
