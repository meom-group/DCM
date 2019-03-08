# DCM User Manual

## DCM fundamentals :
DCM is made of two parts : One for the code compilation (DCMTOOLS) and another one for the code execution (RUNTOOLS).
They are developped side by side and it make sense to have them in a single git repository, under 2 different directories.
In order to ease the synchronisation of the two parts, a well defined hierarchy of directories and files is used throughout the development.
Tools (mainly bash scripts) are provided in order to performs actions in DCM. All these tools have a name starting with dcm\_ and when called 
without argument, a clear *USAGE* message is displayed on the screen.

## [Getting Started](./dcm_getting_started.md)

## [Using DCM for compilation and code maintenance](./dcm_compil_manual.md)

## [Using DCM at runtime](./dcm_rt_manual.md)

## [DCM toolkit](./dcm_toolkit.md) helping/monitoring run production.

