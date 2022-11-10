Tutorial (Docker)
=================

Please download the files in https://github.com/rlrq/MINORg/tree/master/examples. In all the examples below, you should replace "/path/to" with the appropriate full path name, which is usually to the directory containing these example files.

Running MINORg image as container
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the MINORg image as a container using:

.. code-block::
   
   docker run -ti rlrq/minorg

This starts a linux bash session in a directory called '/minorg_docker', so you should see a command prompt like this::

  root@1b2a51c79131:/minorg_docker#

The cryptic alphanumeric string between @ and : is the container ID (in this example, it's ``1b2a51c79131``). Docker will also automatically assign a more human-readable name to the container. To find this human-readable name, start a new terminal session (not in a Docker container) and execute:

.. code-block::

   docker ps

You should see something like this::

  CONTAINER ID   IMAGE     COMMAND   CREATED          STATUS          PORTS     NAMES
  1b2a51c79131   minorg    "bash"    31 seconds ago   Up 29 seconds             great_hodgkin

The container name can be found in the last column.

Alternatively, you can also give container a name by running the image using:

.. code-block::
   
   docker run -ti --name my_minorg rlrq/minorg

In the above example, the newly created container's name will be ``my_minorg``.


Importing files
~~~~~~~~~~~~~~~

As a Docker container is isolated from the host system, files need to be imported into the container, which can be achieved using ``docker cp``. To copy a single file, start a new terminal session (not in a Docker container) and execute:

.. code-block::
   
   docker cp /path/to/file.txt <container ID or name>:/my_minorg/path/to/destination.txt

To copy a directory, use:

.. code-block::
   
   docker cp /path/to/directory <container ID or name>:/my_minorg/path/to/destination


Running MINORg
~~~~~~~~~~~~~~

At this point, you can either follow the command line tutorial starting from :ref:`Tutorial_cli:Setting up the tutorial` to run MINORg in the command line of the container, or start an interactive Python session (Python 3.9 is installed in the container) and follow the Python tutorial starting from :ref:`Tutorial_py:Setting up the tutorial` to run MINORg in Python.


Exporting output files
~~~~~~~~~~~~~~~~~~~~~~

As previously mentioned, the Docker container is isolated from the host system. You may export the files output by MINORg in the container to your system using ``docker cp`` again. To export a single file, use:

.. code-block::
   
   docker cp <container ID or name>:/my_minorg/path/to/file.txt /path/to/destination.txt

To export a directory, use:

.. code-block::
   
   docker cp <container ID or name>:/my_minorg/path/to/directory /path/to/destination


Exiting the container
~~~~~~~~~~~~~~~~~~~~~

When you're done, you may exit the Docker container using ctrl-D or by typing ``exit``.


Reusing the container
~~~~~~~~~~~~~~~~~~~~~

Assuming that you did not remove the container, you can restart an exited Docker container session using:

.. code-block::
   
   docker start -i <container ID or name>

Your files will still be there.


Removing the container
~~~~~~~~~~~~~~~~~~~~~~

If you do not wish to reuse the container, you can delete it using:

.. code-block::

   docker rm <container ID or name>


Docker docs
~~~~~~~~~~~

For more on Docker, you may check out Docker's website: https://docs.docker.com/
