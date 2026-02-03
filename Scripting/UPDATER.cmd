.\..\Library\PortableGit\bin\git.exe clone https://github.com/mkonigs/emma_toolbox.git .\..\Library\GitHub-Repo\

robocopy ..\Library\GitHub-Repo\Scripting ..\Scripting\ /s
robocopy ..\Library\GitHub-Repo\Tests ..\Tests\ /s

timeout 20