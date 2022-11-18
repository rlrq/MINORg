#!/bin/bash

DIR=/mnt/chaelab/rachelle/scripts/minorgpy

cd ${DIR}

## build wheel
BUILD_PYTHON=0
## push latest wheel to (test)pypi
PUSH_PYTHON=0
## update Dockerfile with latest wheel
UPDATE_DOCKERFILE=0
## build docker image
BUILD_DOCKER=0
## push latest docker image to Docker Hub
PUSH_DOCKER=0

while (( "$#" )); do
    case "$1" in
        ## full build + update
        --all|--full) BUILD_PYTHON=1; PUSH_PYTHON=1;
                      UPDATE_DOCKERFILE=1; BUILD_DOCKER=1; PUSH_DOCKER=1;;
        ## standard options
        --build-python) BUILD_PYTHON=1;;
        --push-python) PUSH_PYTHON=1;;
        --update-dockerfile) UPDATE_DOCKERFILE=1;;
        --build-docker) BUILD_DOCKER=1;;
        --push-docker) PUSH_DOCKER=1;;
        ## combined options (platform)
        --python) BUILD_PYTHON=1; PUSH_PYTHON=1;;
        --docker) UPDATE_DOCKERFILE=1; BUILD_DOCKER=1; PUSH_DOCKER=1;;
        ## combined options (actions)
        --build) BUILD_PYTHON=1; UPDATE_DOCKERFILE=1; BUILD_DOCKER=1;;
        --push) PUSH_PYTHON=1; PUSH_DOCKER=1;;
        # -h|--help) man -l ${SCRIPT_DIR}/MANUAL_get_seqs.1; exit 0;;
        # --readme) cat ${SCRIPT_DIR}/README; exit 0;;
    esac
    shift
done

# execute=( "update-dockerfile" )
# # execute=( "python-build" "python-push" "update-dockerfile" "docker-build" "docker-push" )

# execute_options=(
#     "python-build" ## build wheel
#     "python-push" ## push latest wheel to (test)pypi
#     "update-dockerfile" ## update Dockerfile with latest wheel
#     "docker-build" ## build docker image
#     "docker-push" ## push latest docker image to docker hub
# )

## build wheel
if [[ ${BUILD_PYTHON} == 1 ]]; then
    echo "[1] Building wheel"
    python3 -m build
fi

## get name of wheel
echo "[2] Getting name of wheel"
whl=$(ls -t ${DIR}/dist/*.whl | head -1)
whl_basename=$(basename ${whl})
whl_version=$(grep -Po '(?<=minorg-)[^-]+' <<< ${whl_basename})
echo "--wheel is: ${whl_basename}"

## update (test)pypi
if [[ ${PUSH_PYTHON} == 1 ]]; then
    echo "[3] Uploading wheel and source to (test)pypi"
    distr_prefix=${whl_basename%*-py3-*.whl}
    python3 -m twine upload --repository testpypi dist/${distr_prefix}*
fi

## update Dockerfile with name of latest wheel
if [[ ${UPDATE_DOCKERFILE} == 1 ]]; then
    echo "[4] Updating Dockerfile"
    ${DIR}/update_dockerfile.py ${whl} ${whl_version}
fi

## regenerate docker image
if [[ ${BUILD_DOCKER} == 1 ]]; then
    echo "[5.1] Building docker image (full)"
    cp Dockerfile-full Dockerfile
    docker build -t minorg:latest -t minorg:${whl_version} ${DIR}
    echo "[5.1] Building docker image (lite)"
    cp Dockerfile-lite Dockerfile
    docker build -t minorg-lite:latest -t minorg-lite:${whl_version} ${DIR}
fi

## upload wheel to docker hub
if [[ ${PUSH_DOCKER} == 1 ]]; then
    echo "[6] Pushing image to Docker Hub"
    echo "password requirement: >=9 char"
    docker login -u rlrq
    echo "Pushing (full)"
    docker tag minorg:${whl_version} rlrq/minorg:${whl_version}
    docker tag minorg:${whl_version} rlrq/minorg:latest
    docker push rlrq/minorg:${whl_version}
    docker push rlrq/minorg:latest
    echo "Pushing (lite)"
    docker tag minorg-lite:${whl_version} rlrq/minorg-lite:${whl_version}
    docker tag minorg-lite:${whl_version} rlrq/minorg-lite:latest
    docker push rlrq/minorg-lite:${whl_version}
    docker push rlrq/minorg-lite:latest
    docker logout
fi
