package com.lucidlogic.horsetracker.service;

import com.lucidlogic.horsetracker.model.Race;
import com.lucidlogic.horsetracker.model.Racecard;
import com.lucidlogic.horsetracker.model.dto.RaceDTO;
import com.lucidlogic.horsetracker.model.dto.RacecardDTO;

import com.lucidlogic.horsetracker.model.dto.RunnerDTO;
import com.lucidlogic.horsetracker.utils.BeanTransformers;
import com.parse.ParseException;
import com.parse.ParseQuery;

import org.junit.Assert;
import org.junit.Test;

import java.util.Date;
import java.util.List;


/**
 * Created by ian on 10/07/16.
 */
public class ParseServiceTest extends ParseTest{

    @Test
    public void getRacecardsTest() throws ParseException {
        ParseQuery parseQuery = ParseQuery.getQuery(RacecardDTO.class);
        List<RacecardDTO> racecardDTOList = parseQuery.find();
        Assert.assertTrue(racecardDTOList.size() > 0);
        RacecardDTO racecardDTO = ParseService.getRacecards(new Date()).toBlocking().first();
        Assert.assertNotNull(racecardDTO);
        Racecard racecard = BeanTransformers.racecardFromRacecardDTO(racecardDTO, false);
        Assert.assertNotNull(racecard);
    }

    @Test
    public void getRaceTest() throws ParseException {

        RaceDTO raceDTO = ParseService.getRace("yhLeEP2v4m").toBlocking().first();
        Assert.assertNotNull(raceDTO);
        Race race = BeanTransformers.raceFromRaceDTO(raceDTO, true);
        Assert.assertNotNull(race);

    }
}
