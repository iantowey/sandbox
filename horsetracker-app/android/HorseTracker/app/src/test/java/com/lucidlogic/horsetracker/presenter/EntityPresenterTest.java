package com.lucidlogic.horsetracker.presenter;

import com.lucidlogic.horsetracker.model.Entity;
import com.lucidlogic.horsetracker.model.dto.EntityDTO;
import com.lucidlogic.horsetracker.presenter.impl.EntityPresenterImpl;
import com.lucidlogic.horsetracker.service.ParseService;
import com.lucidlogic.horsetracker.service.ParseTest;
import com.lucidlogic.horsetracker.utils.BeanTransformers;
import com.lucidlogic.horsetracker.view.EntityView;
import com.parse.ParseException;

import junit.framework.Assert;

import org.junit.Before;
import org.junit.Test;
import org.mockito.Mockito;
import org.powermock.api.mockito.PowerMockito;

import rx.android.schedulers.AndroidSchedulers;
import rx.schedulers.Schedulers;
import timber.log.Timber;

import static org.mockito.Matchers.any;

/**
 * Created by ian on 01/08/16.
 */
//@PrepareForTest({AndroidSchedulers.class})
public class EntityPresenterTest extends ParseTest {

    private EntityPresenter entityPresenter;

    @Before
    public void before(){

        PowerMockito.stub(PowerMockito.method(AndroidSchedulers.class, "mainThread")).toReturn(Schedulers.immediate());
        PowerMockito.stub(PowerMockito.method(Schedulers.class, "io")).toReturn(Schedulers.immediate());

        entityPresenter = new EntityPresenterImpl();
        EntityView raceView = Mockito.mock(EntityView.class);
        Mockito.doAnswer(answer ->{
            Timber.i("View invoked");
            return null;})
                .when(raceView)
                .updateView(any(Entity.class));
        entityPresenter.setView(raceView);
    }

    @Test
    public void updateViewTest() throws ParseException {
        entityPresenter.updateView("2HSiH6uODf");
        EntityDTO entityDTO = ParseService.getEntity("2HSiH6uODf").toBlocking().first();
        Assert.assertNotNull(entityDTO);
        Entity entity = BeanTransformers.entityFromEntityDTO(entityDTO, true);
        Assert.assertNotNull(entity);
    }
}
