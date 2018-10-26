var url = require('url');

module.exports = {
    parse_url: function(url_str){
        return url.parse(url_str, true);
    }
}