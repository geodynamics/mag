This container hosts a built version of MAG.

docker run -it --rm -v $HOME/mag:/home/mag_user/work geodynamics/mag

This command will start the MAG docker image and give you terminal access. Any changes made in the /home/mag_user/work directory will be reflected on the host machine at home/mag.
