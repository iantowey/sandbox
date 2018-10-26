'use strict';

var CronJob = require('cron').CronJob;

module.exports = {

    schedule_task: function(cron_exp, task, done){
        var job = new CronJob(
            {
                cronTime: cron_exp,
                onTick: task(),
                onComplete: done(),
                start: false,
                timeZone: 'Europe/London'
            });
        return job;
    }
}

// var job = new CronJob(
//     {
//         cronTime: '*/2 * * * * *',
//         onTick: function () {
//             console.log("Done some processing " + new Date().toString());
//         },
//         onComplete: function(){
//             console.log("Job Complete " + new Date().toString());
//         }
//         start: false,
//         timeZone: 'America/Los_Angeles'
//     });
//
// job.start();
// setTimeout(function () {
//     job.stop();
//     console.log('job stopped');
// }, 10000);
//
