# dropout function definition
function try {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
        echo "error with $1" >&2
		exit $status
    fi
    return $status
}

LIB_NAME=libPrimitect
INSTALL_DIR=$LIB_NAME-install
SRC_DIR=.

TMP_ROOT_DIR=$(pwd)
export Eigen3_DIR=$TMP_ROOT_DIR/thirdparty/eigen-install-ios
export Ceres_DIR=$TMP_ROOT_DIR/thirdparty/ceres-install-ios

rm -rf $INSTALL_DIR
mkdir $INSTALL_DIR

try cd $TMP_ROOT_DIR
try ${POLLY_ROOT}/bin/build.py --home $SRC_DIR --toolchain=$IOSTOOLCHAIN --ios-multiarch --config Release --target Primitect
#try ${POLLY_ROOT}/bin/build.py --home $LIB_NAME --toolchain=$IOSTOOLCHAIN --ios-multiarch --config Debug

#try cp -R $LIB_NAME/lib $INSTALL_DIR
#cp -R $LIB_NAME/bin $INSTALL_DIR

try echo "{ \"version\":\"${DATEBLD}\", \"commit\":\"${COMMITSHA}\" }" >> $INSTALL_DIR/Version_primitect.txt


try tar -zcvf $INSTALL_DIR-ios.tar.gz $INSTALL_DIR

# allow upload of latest build to fail
# curl -T $INSTALL_DIR-ios.tar.gz  ftp://slamattack.icg.tugraz.at/latest/ --user ftpuser:ftpuser

INSTALL_DIR=$LIB_NAME-install
FINALLIBNAME=$LIB_NAME.framework/$LIB_NAME
mkdir -p $INSTALL_DIR/lib/$LIB_NAME.framework
LIPSTR="lipo -output $LIB_NAME-install/lib/$FINALLIBNAME -create _builds/$IOSTOOLCHAIN/detect_objects/Release-iphoneos/$LIB_NAME.dylib"
`$LIPSTR`
install_name_tool -id @rpath/$FINALLIBNAME $LIB_NAME-install/lib/$FINALLIBNAME
cp Info.plist $INSTALL_DIR/lib/$LIB_NAME.framework/
exit 0
