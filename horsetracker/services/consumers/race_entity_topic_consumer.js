var kafka = require('kafka-node'),
    HighLevelConsumer = kafka.HighLevelConsumer,
    client = new kafka.Client(),
    consumer = new HighLevelConsumer(
        client,
        [
            { topic: 'topic2' }
        ],
        {
            groupId: 'my-group'
        }
    );

consumer.setOffset('topic2', 0, 0);

consumer.on('message', function (message) {
    console.log('got message ' + JSON.stringify(message));
});