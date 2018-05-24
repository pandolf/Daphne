voms-proxy-init --voms cms --valid 168:00
cp /tmp/x509up_u27697 `pwd`/X509_USER_PROXY
export X509_USER_PROXY=`pwd`/X509_USER_PROXY
echo "X509_USER_PROXY: "
echo $X509_USER_PROXY
