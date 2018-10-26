'use strict';

var dateFormat = require('dateformat');


module.exports = {
    sp_raceards_url: function(config,dt) {
        return config.sporting_life_com.url 
            + config.sporting_life_com.racecard.url.replace('{dt}', dateFormat(dt, config.sporting_life_com.racecard.url_dt_path_param_format));
    }
};