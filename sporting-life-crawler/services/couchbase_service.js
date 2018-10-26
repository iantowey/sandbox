'use strict';

var couchbase = require("couchbase");
var conf = require('../config/' + process.env.NODE_ENV + '.json');

var cluster = new couchbase.Cluster(conf.couchbase.server);
var bucket = cluster.openBucket(conf.couchbase.bucket, conf.couchbase.bucket_password, function(err) {
    if (err) {
        console.error("Fatal Error : cant connect to couchbase");
        process.exit();
    }
});

bucket.operationTimeout = conf.couchbase.connection_timeout;

module.exports = {
    get_bucket : function(){
        return bucket;
    }
}