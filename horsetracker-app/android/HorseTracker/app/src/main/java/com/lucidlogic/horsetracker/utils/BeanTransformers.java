package com.lucidlogic.horsetracker.utils;


import com.annimon.stream.Collectors;
import com.annimon.stream.Stream;
import com.lucidlogic.horsetracker.model.Entity;
import com.lucidlogic.horsetracker.model.EntityDetailFuture;
import com.lucidlogic.horsetracker.model.EntityDetailHistorical;
import com.lucidlogic.horsetracker.model.Meeting;
import com.lucidlogic.horsetracker.model.Race;
import com.lucidlogic.horsetracker.model.Racecard;
import com.lucidlogic.horsetracker.model.Runner;
import com.lucidlogic.horsetracker.model.dto.EntityDTO;
import com.lucidlogic.horsetracker.model.dto.EntityDetailFutureDTO;
import com.lucidlogic.horsetracker.model.dto.EntityDetailHistoricalDTO;
import com.lucidlogic.horsetracker.model.dto.MeetingDTO;
import com.lucidlogic.horsetracker.model.dto.RaceDTO;
import com.lucidlogic.horsetracker.model.dto.RacecardDTO;
import com.lucidlogic.horsetracker.model.dto.RunnerDTO;

/**
 * Created by ian on 08/06/16.
 */
public class BeanTransformers {

    public static Racecard racecardFromRacecardDTO(RacecardDTO racecardDTO, boolean includeRunners) {
        Racecard racecard = new Racecard();
        racecard.setId(racecardDTO.getObjectId());
        racecard.getMeetings().addAll(
                Stream.of(racecardDTO.getMeetings())
                        .map(meetingDTO -> meetingFromMeetingDTO(meetingDTO, includeRunners))
                        .collect(Collectors.toList())
        );
        racecard.setDate(racecardDTO.getDate());


        return racecard;
    }

    public static Meeting meetingFromMeetingDTO(MeetingDTO meetingDTO, boolean includeRunners) {
        Meeting meeting = new Meeting();

        meeting.setId(meetingDTO.getObjectId());
        meeting.setGoing(meetingDTO.getTrackGoing());
        meeting.setSurface(meetingDTO.getTrackSurface());
        meeting.setTrack(meetingDTO.getTrack());

        meeting.getRaces().addAll(
                Stream.of(meetingDTO.getRaces())
                .map(raceDTO -> raceFromRaceDTO(raceDTO, includeRunners))
                .collect(Collectors.toList()));
        return meeting;
    }

    public static Race raceFromRaceDTO(RaceDTO raceDTO, boolean includeRunners) {
        Race race = new Race();
        race.setId(raceDTO.getObjectId().trim());
        race.setAbandoned(raceDTO.getAbandoned().trim());
        race.setName(raceDTO.getName().trim());
        race.setTime(raceDTO.getTime().trim());
        if(includeRunners && raceDTO.getRunners() != null && raceDTO.getRunners().size() > 0 ) {
            race.getRunners().addAll(
                    Stream.of(raceDTO.getRunners())
                            .map(runnerDTO -> runnerFromRunnerDTO(runnerDTO))
                            .collect(Collectors.toList())
            );
        }
        return race;
    }

    public static Runner runnerFromRunnerDTO(RunnerDTO runnerDTO) {
        Runner runner = new Runner ();
        runner.setId(runnerDTO.getObjectId().trim());
        runner.setAge(runnerDTO.getAge().trim());
        runner.setBreeding(runnerDTO.getBreeding().trim());
        runner.setCreatedAt(runnerDTO.getCreatedAt());
        runner.setForm(runnerDTO.getForm());
        runner.setOr(runnerDTO.getOr());
        runner.setRunning(runnerDTO.isRunning());
        runner.setSilkImgLink(runnerDTO.getSilkImgLink());
        runner.setStall(runnerDTO.getStall());
        runner.setHorse(entityFromEntityDTO(runnerDTO.getHorse(), false));
        runner.setTrainer(entityFromEntityDTO(runnerDTO.getTrainer(), false));
        runner.setJockey(entityFromEntityDTO(runnerDTO.getJockey(), false));
        return runner;
    }

    public static Entity entityFromEntityDTO(EntityDTO entityDTO, boolean details) {
        Entity entity = new Entity();
        entity.setName(entityDTO.getName());
        entity.setId(entityDTO.getObjectId());
        entity.setProfileUrl(entityDTO.getProfileUrl());
        entity.setType(entityDTO.getType());
        if(details) {
            entity.setEntityDetailFuture(entityDetailFutureFromEntityDetailFutureDTO(entityDTO.getEntityDetailFutureDTO()));
            entity.setEntityDetailHistorical(entityDetailHistoricalFromEntityDetailHistoricalDTO(entityDTO.getEntityDetailHistoricalDTO()));
        }
        return entity;
    }

    private static EntityDetailHistorical entityDetailHistoricalFromEntityDetailHistoricalDTO(EntityDetailHistoricalDTO entityDetailHistoricalDTO) {

        EntityDetailHistorical entityDetailHistorical = new EntityDetailHistorical();
        entityDetailHistorical.setId(entityDetailHistoricalDTO.getObjectId());
        entityDetailHistorical.setName(entityDetailHistoricalDTO.getName());
        entityDetailHistorical.setAge(entityDetailHistoricalDTO.getAge());
        entityDetailHistorical.setCollectionDate(entityDetailHistoricalDTO.getCollectionDate());
        entityDetailHistorical.setDam(entityDetailHistoricalDTO.getDam());
        entityDetailHistorical.setSire(entityDetailHistoricalDTO.getSire());
        entityDetailHistorical.setHistoricalForm(entityDetailHistoricalDTO.getHistoricalForm());
        entityDetailHistorical.setHistoricalSummary(entityDetailHistoricalDTO.getHistoricalSummary());
        entityDetailHistorical.setOwner(entityDetailHistoricalDTO.getOwner());

        return entityDetailHistorical;
    }

    private static EntityDetailFuture entityDetailFutureFromEntityDetailFutureDTO(EntityDetailFutureDTO entityDetailFutureDTO) {
        EntityDetailFuture entityDetailFuture = new EntityDetailFuture();

        entityDetailFuture.setId(entityDetailFutureDTO.getObjectId());
        entityDetailFuture.setName(entityDetailFutureDTO.getName());
        entityDetailFuture.setAge(entityDetailFutureDTO.getAge());
        entityDetailFuture.setCollectionDate(entityDetailFutureDTO.getCollectionDate());
        entityDetailFuture.setDam(entityDetailFutureDTO.getDam());
        entityDetailFuture.setSire(entityDetailFutureDTO.getSire());
        entityDetailFuture.setOwner(entityDetailFutureDTO.getOwner());
        entityDetailFuture.setFutEnt(entityDetailFutureDTO.getFutEnt());

        return entityDetailFuture;
    }
}
